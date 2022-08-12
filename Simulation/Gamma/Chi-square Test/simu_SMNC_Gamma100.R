#!/bin/bash
#SBATCH --job-name=SMNCG100
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=25
#SBATCH --cluster=mpi
#SBATCH --partition=opa-high-mem
#SBATCH --mail-user=zhw61@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00

module purge
module load gcc/8.2.0
module load r/3.6.0

#This is a comment

export OMPI_MCA_rmaps_base_oversubscribe=1

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl="self,tcp"

#cp $SLURM_SUBMIT_DIR/* $SLURM_SCRATCH
#cd $SLURM_SCRATCH

cd $SLURM_SUBMIT_DIR

#echo $SLURM_JOB_NODELIST > slurm_nodefile.txt

scontrol show hostnames | awk "{for(i=0;i<$SLURM_NTASKS_PER_NODE;i++)print}" > node_list.txt

# Telling cluster that you are using R
R --vanilla > mysnow.out <<EOF

#################################################################################
# Loading snowfall and looking for available nodes, then initializing:
#################################################################################
Sys.setenv(OMPI_MCA_btl="tcp,self")

library(snowfall)
##pbsnodefile = Sys.getenv("PBS_NODEFILE")
##pbsnodefile = Sys.getenv("SLURM_JOB_NODELIST")
##machines <- scan(pbsnodefile, what="")

machines <- scan("$SLURM_SUBMIT_DIR/node_list.txt", what="")
machines
nmach = length(machines)

sfInit(parallel=TRUE,type='MPI',cpus=nmach,socketHosts=machines)

#################################################################################
# Loading Other necessary R Libraries:
#################################################################################
sfLibrary(survival)
sfLibrary(lme4)

#################################################################################
#  Load External Data
#################################################################################


#################################################################################
#  Main function
#################################################################################



#################################################################################
#  Exporting the data and functions needed by the worker nodes:
#################################################################################
sfExportAll()
#################################################################################
#  Creating a wrapper function
#################################################################################


wrapper <- function(name){
  melt=function(data,id){
    trans=names(data)[!names(data)%in%id]
    temp=data[rep(1:nrow(data),length(trans)),id]
    temp[,'value']=c(as.matrix(data[,trans]))
    temp[,'variable']=rep(trans,each=nrow(data))
    return(temp)
  }
  nvisitf=function(lambda,TT){
    t=TT
    UU=VV=NULL
    nvisit=NULL
    k=p=0
    for(i in 1:length(lambda))
      while (TRUE) {
        temp=ifelse(p==0,0,rexp(1,lambda[i]))
        if(k+temp>t[i]){
          nvisit=c(nvisit,p)
          VV=c(VV,k)
          k=p=0
          break
        } else {
          k=k+temp
          UU=c(UU,k)
          p=p+1
        }
      }
    OO=NULL
    OO[['time']]=UU
    OO[['nvisit']]=nvisit
    OO[['last']]=VV
    return(OO)
  }
  ndim=1000
  beta=c(0.2,0.2,0.2,0.2,-0.2)
  z1=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z2=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z3=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z4=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z5=runif(ndim)
  beta0=3
  sigma=.1
  Z=cbind(z1,z2,z3,z4,z5)
  U=runif(ndim)
  w = log(-log(1-U))					# inverse cdf of extreme value dist
  TT = exp(beta0 + Z%*%beta + sigma*w)		# Weibull model
  X=runif(ndim,0,20)
  St=pmin(TT,X)
  delta=ifelse(TT<X,1,0)
  
  z12=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z22=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z32=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z42=pmin(pmax(rnorm(ndim),-3.5),3.5)
  z52=runif(ndim)
  Z2=cbind(z12,z22,z32,z42,z52)
  U2=runif(ndim)
  w2 = log(-log(1-U2))					# inverse cdf of extreme value dist
  TT2 = exp(beta0 + Z2%*%beta + sigma*w2)		# Weibull model
  X2=runif(ndim,0,20)
  St2=pmin(TT2,20-X2)
  delta2=ifelse(TT2<20-X2,1,0)
  
  
  gamma=c(rep(0.1,4),-0.1)
  gamma0=0
  lambda=exp(gamma0 + Z%*%gamma)
  group1=nvisitf(lambda,St)
  visit.num=group1[['nvisit']]
  lambda2=exp(gamma0 + Z2%*%gamma)
  group2=nvisitf(lambda2,St2)
  visit.num2=group2[['nvisit']]
  
  library(MASS)
  
  
  NP=NP2=NULL
  for(i in 1:ndim){
    np.means=rep(0,6*visit.num[i])
    np.covs=matrix(rep(0,6*6*visit.num[i]^2),ncol=6*visit.num[i])
    for(k in 1:(6*visit.num[i])){
      for(j in 1:(6*visit.num[i])){
        if(k==j)
          np.covs[k,j]=100
        else if((k-j)%%6==0)
          np.covs[k,j]=60
        else if((k-1)%/%6==(j-1)%/%6)
          np.covs[k,j]=25
        else
          np.covs[k,j]=15
      }
    }
    nps=mvrnorm(1,np.means,np.covs/100)
    nps=150-qgamma(pnorm(nps),100,1)
    for(j in 1:visit.num[i]){
      NP=rbind(NP,c(j,i,nps[((j-1)*6+1):(j*6)]))
    }
  }
  NP=as.data.frame(NP)
  names(NP)=c('visit','subid','np1','np2','np3','np4','np5','np6')
  NP[,'t']=group1[['time']]
  NP[,'t2']=NP[,'t']^2/10
  NP[,'t3']=NP[,'t']^3/100
  
  survdata=as.data.frame(cbind(z1,z2,z3,z4,z5,St,delta))
  survdata2=as.data.frame(cbind(z12,z22,z32,z42,z52,St2,delta2))
  names(survdata)=c('z1','z2','z3','z4','z5','St','event')
  names(survdata2)=c('z1','z2','z3','z4','z5','St','event')
  cox = coxph(Surv(St, delta) ~ z1+z2+z3+z4+z5,data=survdata)
  exp_survtime=summary(survfit(cox,survdata2))[['table']][,'median']
  exp_survtime[is.na(exp_survtime)]=20
  
  countdata=as.data.frame(cbind(z1,z2,z3,z4,z5,St,visit.num-1))
  countdata2=as.data.frame(cbind(z12,z22,z32,z42,z52,St2,visit.num2-1))
  names(countdata)=c('z1','z2','z3','z4','z5','St','visits')
  names(countdata2)=c('z1','z2','z3','z4','z5','St','visits')
  count <- glm(visits~z1+z2+z3+z4+z5+offset(log(St)),family=poisson(link=log),data=countdata)
  lambdas=exp(as.matrix(cbind(rep(1,ndim),countdata2[,1:5]))%*%as.matrix(count[['coefficients']])[,1])
  
  KK=melt(NP,id=c('visit','subid','t','t2','t3'))
  mlme=summary(lmer(value~factor(variable)+factor(variable)*t+factor(variable)*t2+factor(variable)*t3+(1|subid)+(1|visit:subid)+(1|variable:subid),data=KK))
  
  NP.linear=NP
  NP.linear[,c('np1','np2','np3')]=NP.linear[,c('np1','np2','np3')]-0.6*NP.linear[,'t']
  NP.linear[,c('np4','np5','np6')]=NP.linear[,c('np4','np5','np6')]-0.8*NP.linear[,'t']
  KK=melt(NP.linear,id=c('visit','subid','t','t2','t3'))
  mlme.linear=summary(lmer(value~factor(variable)+factor(variable)*t+factor(variable)*t2+factor(variable)*t3+(1|subid)+(1|visit:subid)+(1|variable:subid),data=KK))
  
  NP.quad=NP
  NP.quad[,c('np1','np2','np3')]=NP.quad[,c('np1','np2','np3')]-0.08*NP.quad[,'t']^2+0.2*NP.quad[,'t']
  NP.quad[,c('np4','np5','np6')]=NP.quad[,c('np4','np5','np6')]-0.06*NP.quad[,'t']^2+0.1*NP.quad[,'t']
  KK=melt(NP.quad,id=c('visit','subid','t','t2','t3'))
  mlme.quad=summary(lmer(value~factor(variable)+factor(variable)*t+factor(variable)*t2+factor(variable)*t3+(1|subid)+(1|visit:subid)+(1|variable:subid),data=KK))
  
  NP.cubic=NP
  NP.cubic[,c('np1','np2','np3')]=NP.cubic[,c('np1','np2','np3')]-0.008*NP.cubic[,'t']^3+0.08*NP.cubic[,'t']^2+0.55*NP.cubic[,'t']
  NP.cubic[,c('np4','np5','np6')]=NP.cubic[,c('np4','np5','np6')]-0.007*NP.cubic[,'t']^3+0.06*NP.cubic[,'t']^2+0.55*NP.cubic[,'t']
  KK=melt(NP.cubic,id=c('visit','subid','t','t2','t3'))
  mlme.cubic=summary(lmer(value~factor(variable)+factor(variable)*t+factor(variable)*t2+factor(variable)*t3+(1|subid)+(1|visit:subid)+(1|variable:subid),data=KK))
  
  
  
  estmean=mlme[['coefficients']][1,1]+c(0,mlme[['coefficients']][2:6,1])
  estmean.linear=mlme.linear[['coefficients']][1,1]+c(0,mlme.linear[['coefficients']][2:6,1])
  estmean.quad=mlme.quad[['coefficients']][1,1]+c(0,mlme.quad[['coefficients']][2:6,1])
  estmean.cubic=mlme.cubic[['coefficients']][1,1]+c(0,mlme.cubic[['coefficients']][2:6,1])
  
  cont1=mlme[['coefficients']][7,1]+c(0,mlme[['coefficients']][10:14,1])
  cont2=mlme[['coefficients']][8,1]+c(0,mlme[['coefficients']][15:19,1])
  cont3=mlme[['coefficients']][9,1]+c(0,mlme[['coefficients']][20:24,1])
  
  linear1=mlme.linear[['coefficients']][7,1]+c(0,mlme.linear[['coefficients']][10:14,1])
  linear2=mlme.linear[['coefficients']][8,1]+c(0,mlme.linear[['coefficients']][15:19,1])
  linear3=mlme.linear[['coefficients']][9,1]+c(0,mlme.linear[['coefficients']][20:24,1])
  
  quad1=mlme.quad[['coefficients']][7,1]+c(0,mlme.quad[['coefficients']][10:14,1])
  quad2=mlme.quad[['coefficients']][8,1]+c(0,mlme.quad[['coefficients']][15:19,1])
  quad3=mlme.quad[['coefficients']][9,1]+c(0,mlme.quad[['coefficients']][20:24,1])
  
  cubic1=mlme.cubic[['coefficients']][7,1]+c(0,mlme.cubic[['coefficients']][10:14,1])
  cubic2=mlme.cubic[['coefficients']][8,1]+c(0,mlme.cubic[['coefficients']][15:19,1])
  cubic3=mlme.cubic[['coefficients']][9,1]+c(0,mlme.cubic[['coefficients']][20:24,1])
  
  
  res=NULL
  for(i in 1:ndim){
    time=group2[['time']][(sum(c(0,group2[['nvisit']])[1:i])+1):sum(group2[['nvisit']][1:i])]
    np.means2=rep(0,6*visit.num2[i])
    np.covs2=matrix(rep(0,6*6*visit.num2[i]^2),ncol=6*visit.num2[i])
    for(k in 1:(6*visit.num2[i])){
      for(j in 1:(6*visit.num2[i])){
        if(k==j)
          np.covs2[k,j]=100
        else if((k-j)%%6==0)
          np.covs2[k,j]=60
        else if((k-1)%/%6==(j-1)%/%6)
          np.covs2[k,j]=25
        else
          np.covs2[k,j]=15
      }
    }
    nps2=mvrnorm(1,np.means2,np.covs2/100)
    nps2=150-qgamma(pnorm(nps2),100,1)
    nps2.linear=nps2-rep(c(0.6,0.6,0.6,0.8,0.8,0.8),visit.num2[i])*rep(time,each=6)
    nps2.quad=nps2+rep(c(0.2,0.2,0.2,0.1,0.1,0.1),visit.num2[i])*rep(time,each=6)-rep(c(0.8,0.8,0.8,0.6,0.6,0.6),visit.num2[i])*rep(time^2/10,each=6)
    nps2.cubic=nps2+rep(c(0.55,0.55,0.55,0.55,0.55,0.55),visit.num2[i])*rep(time,each=6)+rep(c(0.8,0.8,0.8,0.6,0.6,0.6),visit.num2[i])*rep(time^2/10,each=6)-rep(c(0.8,0.8,0.8,0.7,0.7,0.7),visit.num2[i])*rep(time^3/100,each=6)
    
    MEANS=rep(estmean,visit.num2[i])+rep(cont1,visit.num2[i])*rep(time,each=6)+rep(cont2,visit.num2[i])*rep(time^2/10,each=6)+rep(cont3,visit.num2[i])*rep(time^3/100,each=6)
    estcov=matrix(rep(0,6*6*visit.num2[i]^2),ncol=6*visit.num2[i])
    for(k in 1:(6*visit.num2[i])){
      for(j in 1:(6*visit.num2[i])){
        if(k==j)
          estcov[k,j]=sum(as.numeric(mlme[['varcor']]))+mlme[['sigma']]^2
        else if((k-j)%%6==0)
          estcov[k,j]=sum(as.numeric(mlme[['varcor']])[-1])
        else if((k-1)%/%6==(j-1)%/%6)
          estcov[k,j]=sum(as.numeric(mlme[['varcor']])[-2])
        else
          estcov[k,j]=sum(as.numeric(mlme[['varcor']])[3])
      }
    }
    
    MEANS.linear=rep(estmean.linear,visit.num2[i])+rep(linear1,visit.num2[i])*rep(time,each=6)+rep(linear2,visit.num2[i])*rep(time^2/10,each=6)+rep(linear3,visit.num2[i])*rep(time^3/100,each=6)
    estcov.linear=matrix(rep(0,6*6*visit.num2[i]^2),ncol=6*visit.num2[i])
    for(k in 1:(6*visit.num2[i])){
      for(j in 1:(6*visit.num2[i])){
        if(k==j)
          estcov.linear[k,j]=sum(as.numeric(mlme.linear[['varcor']]))+mlme.linear[['sigma']]^2
        else if((k-j)%%6==0)
          estcov.linear[k,j]=sum(as.numeric(mlme.linear[['varcor']])[-1])
        else if((k-1)%/%6==(j-1)%/%6)
          estcov.linear[k,j]=sum(as.numeric(mlme.linear[['varcor']])[-2])
        else
          estcov.linear[k,j]=sum(as.numeric(mlme.linear[['varcor']])[3])
      }
    }
    
    MEANS.quad=rep(estmean.quad,visit.num2[i])+rep(quad1,visit.num2[i])*rep(time,each=6)+rep(quad2,visit.num2[i])*rep(time^2/10,each=6)+rep(quad3,visit.num2[i])*rep(time^3/100,each=6)
    estcov.quad=matrix(rep(0,6*6*visit.num2[i]^2),ncol=6*visit.num2[i])
    for(k in 1:(6*visit.num2[i])){
      for(j in 1:(6*visit.num2[i])){
        if(k==j)
          estcov.quad[k,j]=sum(as.numeric(mlme.quad[['varcor']]))+mlme.quad[['sigma']]^2
        else if((k-j)%%6==0)
          estcov.quad[k,j]=sum(as.numeric(mlme.quad[['varcor']])[-1])
        else if((k-1)%/%6==(j-1)%/%6)
          estcov.quad[k,j]=sum(as.numeric(mlme.quad[['varcor']])[-2])
        else
          estcov.quad[k,j]=sum(as.numeric(mlme.quad[['varcor']])[3])
      }
    }
    
    MEANS.cubic=rep(estmean.cubic,visit.num2[i])+rep(cubic1,visit.num2[i])*rep(time,each=6)+rep(cubic2,visit.num2[i])*rep(time^2/10,each=6)+rep(cubic3,visit.num2[i])*rep(time^3/100,each=6)
    estcov.cubic=matrix(rep(0,6*6*visit.num2[i]^2),ncol=6*visit.num2[i])
    for(k in 1:(6*visit.num2[i])){
      for(j in 1:(6*visit.num2[i])){
        if(k==j)
          estcov.cubic[k,j]=sum(as.numeric(mlme.cubic[['varcor']]))+mlme.cubic[['sigma']]^2
        else if((k-j)%%6==0)
          estcov.cubic[k,j]=sum(as.numeric(mlme.cubic[['varcor']])[-1])
        else if((k-1)%/%6==(j-1)%/%6)
          estcov.cubic[k,j]=sum(as.numeric(mlme.cubic[['varcor']])[-2])
        else
          estcov.cubic[k,j]=sum(as.numeric(mlme.cubic[['varcor']])[3])
      }
    }
    
    LMNCs=NULL
    directions=NULL
    for(k in 1:visit.num2[i]){
      directions=rbind(directions,c(ifelse(sum(nps2[(k*6-5):(k*6)]-MEANS[(k*6-5):(k*6)])<0,-1,1),
                                    ifelse(sum(nps2.linear[(k*6-5):(k*6)]-MEANS.linear[(k*6-5):(k*6)])<0,-1,1),
                                    ifelse(sum(nps2.quad[(k*6-5):(k*6)]-MEANS.quad[(k*6-5):(k*6)])<0,-1,1),
                                    ifelse(sum(nps2.cubic[(k*6-5):(k*6)]-MEANS.cubic[(k*6-5):(k*6)])<0,-1,1)))
      LMNCs=rbind(LMNCs,c((nps2[1:(k*6)]-MEANS[1:(k*6)])%*%solve(estcov[1:(k*6),1:(k*6)])%*%(nps2[1:(k*6)]-MEANS[1:(k*6)]),
                          (nps2.linear[1:(k*6)]-MEANS.linear[1:(k*6)])%*%solve(estcov.linear[1:(k*6),1:(k*6)])%*%(nps2.linear[1:(k*6)]-MEANS.linear[1:(k*6)]),
                          (nps2.quad[1:(k*6)]-MEANS.quad[1:(k*6)])%*%solve(estcov.quad[1:(k*6),1:(k*6)])%*%(nps2.quad[1:(k*6)]-MEANS.quad[1:(k*6)]),
                          (nps2.cubic[1:(k*6)]-MEANS.cubic[1:(k*6)])%*%solve(estcov.cubic[1:(k*6),1:(k*6)])%*%(nps2.cubic[1:(k*6)]-MEANS.cubic[1:(k*6)])))
    }
    diffs=LMNCs-rbind(rep(0,4),LMNCs[-visit.num2[i],])
    diffs=diffs*directions
    alphas=seq(0.001,0.1,by=0.0001)
    N_prect=ceiling(min(exp_survtime[i],20-X2[i])*lambdas[i])+1
    
    signal=NULL
    for(alpha in alphas){
      if(visit.num2[i]<=N_prect)
        thres=2*(N_prect-1:visit.num2[i]+1)/N_prect/(N_prect+1)
      else
        thres=c(2*(N_prect:1)/N_prect/(N_prect+1),rep(2/N_prect/(N_prect+1),visit.num2[i]-N_prect))
      
      signal=rbind(signal, ifelse(colSums(-diffs>qchisq(1-2*thres*alpha,6))>0,1,0))
    }
    res=cbind(res,signal)
    cat(i,'\n')
  }
  
  return(cbind(rowMeans(res[,seq(1,by=4,length=ndim)]),rowMeans(res[,seq(2,by=4,length=ndim)]),rowMeans(res[,seq(3,by=4,length=ndim)]),rowMeans(res[,seq(4,by=4,length=ndim)])))
  
}

#################################################################################
#  Running the function nSim times:
#################################################################################
nSim <- 1000
SMNC_Gamma100 <- sfLapply(1:nSim,wrapper)

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(SMNC_Gamma100,file = 'SMNC_Gamma100.rda')

#################################################################################
#  Close all connections:
#################################################################################
sfStop()

EOF

