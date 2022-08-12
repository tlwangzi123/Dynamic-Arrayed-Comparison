#################################################################################
# Loading snowfall and looking for available nodes, then initializing:
#################################################################################
Sys.setenv(OMPI_MCA_btl="tcp,self")

library(snowfall)
##pbsnodefile = Sys.getenv("PBS_NODEFILE")
##pbsnodefile = Sys.getenv("SLURM_JOB_NODELIST")
##machines <- scan(pbsnodefile, what="")

##machines <- scan("$SLURM_SUBMIT_DIR/node_list.txt", what="")
machines <- scan("node_list.txt", what="")
machines
nmach = length(machines)

sfInit(parallel=TRUE,type='MPI',cpus=nmach,socketHosts=machines)
#################################################################################
# Loading Other necessary R Libraries:
#################################################################################
sfLibrary(survival)
sfLibrary(lme4)
sfLibrary(MASS)
sfLibrary(VGAM)
#################################################################################
#  Load External Data
#################################################################################
#functions
#source("/ihome/ycheng/ziw43/Project/Project1/functions/RepNAmin.R")
#data
#load("/ihome/ycheng/baw90/xxx.RData")
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
wrapper <- function(seed){
  tryCatch(
    expr={
    melt=function(data,id){
    trans=names(data)[!names(data)%in%id]
    temp=data[rep(1:nrow(data),length(trans)),id]
    temp[,'value']=c(as.matrix(data[,trans]))
    temp[,'variable']=rep(trans,each=nrow(data))
    return(temp)
  }
  nvisitf_NB=function(delta_0,lambda,TT){
    t=TT
    UU=VV=NULL
    nvisit=NULL
    k=p=0
    for(i in 1:length(lambda)){
    rat = 1
     #rat = rgamma(1, shape = 1/delta_0, rate = 1/delta_0)
      while (TRUE) {
        # rexp(1, lambda[i])
        #temp=ifelse(p==0,0, ((1 - runif(1))^(-delta_0) - 1)/(delta_0*lambda[i]))
        temp = ifelse(p==0, 0, rexp(1, rat * lambda[i]))
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
      }
    OO=NULL
    OO[['time']]=UU
    OO[['nvisit']]=nvisit
    OO[['last']]=VV
    return(OO)
    }
  mu_f = function(r, obj, X){
#  X = Z
  n_obj = dim(obj)[1]
  n_r = length(r)
  r = as.matrix(r)
  X_add = cbind(rep(1, n_obj), X)
  ans = as.numeric(X_add %*% r[1:(n_r-1),])  # the last term of r is alpha
  obj$t * exp(ans)  
}

log_pdf = function(r, mu, obj_num, X){ # not directly calculate pdf since there is a factorial term
  n = length(obj_num)
  n_r = length(r)
  ans = rep(0, n)
#  obj_num = obj[,"num"]
  
  for (i in 1:n){
    for (j in 1:obj_num[i]){
      ans[i] = ans[i] - log(j)
    }
    if (obj_num[i] == 0) ans[i] = 0   # since 0! = 1
    ans[i] = ans[i] + obj_num[i] * log(mu[i] / (1 + r[n_r] * mu[i])) + 
      (obj_num[i] - 1) * log(1 + r[n_r] * obj_num[i]) -
      mu[i] * (1 + r[n_r] * obj_num[i]) / (1 + r[n_r] * mu[i])  
  }
  ans
}

log_surv = function(r, mu, obj){  # includes the boundary y_i
  n_obj = dim(obj)[1]
  ans_surv = rep(1, n_obj)
  
  for (i in 1:n_obj){
    if (obj[i,"num"] >= 1){
      temp = 0
      for (j in 0:(obj[i,"num"] - 1))
        temp = temp + exp(log_pdf(r, mu[i], j))  
      ans_surv[i] = ans_surv[i] - temp     
      # when temp is very small, ans_surv[i] will print the value 1 when called, 
      # although the stored value in ans_surv[i] is not exactly 1
    }
    if (obj[i,"num"] == 0) ans_surv[i] = 1  
  }
  ans_surv = log(ans_surv)
  ans_surv
}

neg_log_lik = function(r, obj, X){ 
  # note: the censoring indicator in the paper is different than that in the survival framework
  # and we use the one in the survival framework (which is 1 - the censoring indicator in the paper)
  n_obj = dim(obj)[1]
  mu = mu_f(r, obj, X)
  log_pdf_obj = log_pdf(r, mu, obj[,"num"], X)
  log_surv_obj = log_surv(r, mu, obj)
  
  ans_lik = 0 
  for (i in 1:n_obj){
    ans_lik = ans_lik + obj[i,"ind"] * log_pdf_obj[i] + 
      (1 - obj[i,"ind"]) * log_surv_obj[i]
  }
  ans_lik = (-1) * ans_lik
  ans_lik
}

ndim=1000
beta=c(0.2,0.2,0.2,0.2,-0.2)
z1=pmin(pmax(rnorm(ndim),-3.5),3.5) # avoid extreme values
z2=pmin(pmax(rnorm(ndim),-3.5),3.5)
z3=pmin(pmax(rnorm(ndim),-3.5),3.5)
z4=pmin(pmax(rnorm(ndim),-3.5),3.5)
z5=runif(ndim)
beta0=3

# Negative-Binomial regression for frequency
delta_0 = 1

sigma=.1
Z=cbind(z1,z2,z3,z4,z5)

# Use log-logistic AFT model for time
U=runif(ndim)
w = log(-log(1-U))
TT = exp(beta0 + Z%*%beta + sigma*w)		# Weibull model
X=runif(ndim,0,20)
St=pmin(TT,X) # enrolled or already dead before enrollment
delta=ifelse(TT<X,1,0)
  
z12=pmin(pmax(rnorm(ndim),-3.5),3.5)
z22=pmin(pmax(rnorm(ndim),-3.5),3.5)
z32=pmin(pmax(rnorm(ndim),-3.5),3.5)
z42=pmin(pmax(rnorm(ndim),-3.5),3.5)
z52=runif(ndim)
Z2=cbind(z12,z22,z32,z42,z52)

U=runif(ndim)
w2 = log(-log(1-U))
TT2 = exp(beta0 + Z2%*%beta + sigma*w2)		# Weibull model
X2=runif(ndim,0,20)
St2=pmin(TT2,20-X2)  
delta2=ifelse(TT2<20-X2,1,0)
# Q: why the data generating mechanism of St1, delat1 and St2, delta2 are not the same?
  
gamma=c(rep(0.1,4),-0.1)
gamma0=0
# eta = rgamma(1, shape = 1/delta, rate = 1/delta)
lambda=exp(gamma0 + Z%*%gamma) # model 5
group1=nvisitf_NB(delta_0, lambda, St)
visit.num=group1[['nvisit']]

#eta2 = rgamma(1, shape = 1/delta, rate = 1/delta)
lambda2=exp(gamma0 + Z2%*%gamma) # model 5
group2=nvisitf_NB(delta_0, lambda2,St2)
visit.num2=group2[['nvisit']]

NP=NP2=NULL
for(i in 1:ndim){
  np.means=rep(50,6*visit.num[i])
  np.covs=matrix(rep(0,6*6*visit.num[i]^2),ncol=6*visit.num[i])
  # construct (true) covariance matrix of cognitive domains  for each subject i
  for(k in 1:(6*visit.num[i])){
    for(j in 1:(6*visit.num[i])){
      if(k==j)
        np.covs[k,j]=100
      else if((k-j)%%6==0)
        np.covs[k,j]=60
      else if((k-1)%/%6==(j-1)%/%6)
        np.covs[k,j]=25
      else
        # cov of different cognitive domains at different visits
        np.covs[k,j]=15
    }
  }
    # generate all 6*n_visit cognitive domain errors for each subject i
    # add constant mean trend
  nps=mvrnorm(1,np.means,np.covs)
    # convert nps into n_visit rows, each row store 6 cognitive domain scores at that visit
    # (from np1 to np6)
    # "i"("subid") column store the id of that subject (i)
    # "j"("visit") column store the visit number of each subject (among the total n_visit visits)
  for(j in 1:visit.num[i]){
    NP=rbind(NP,c(j,i,nps[((j-1)*6+1):(j*6)]))
  }
}
  
NP=as.data.frame(NP)
names(NP)=c('visit','subid','np1','np2','np3','np4','np5','np6')
NP[,'t']=group1[['time']]
NP[,'t2']=NP[,'t']^2/10
NP[,'t3']=NP[,'t']^3/100
  # Q: why not directly use t2=t^2 and t3=t^3?

survdata=as.data.frame(cbind(z1,z2,z3,z4,z5,St,delta))
  survdata2=as.data.frame(cbind(z12,z22,z32,z42,z52,St2,delta2))
  names(survdata)=c('z1','z2','z3','z4','z5','St','event')
  names(survdata2)=c('z1','z2','z3','z4','z5','St','event')
  cox = coxph(Surv(St, delta) ~ z1+z2+z3+z4+z5,data=survdata) # model 4
  # calculate \hat{T_i}
  exp_survtime=summary(survfit(cox,survdata2))[['table']][,'median']
  # 
  exp_survtime[is.na(exp_survtime)]=20
  
  countdata=as.data.frame(cbind(z1,z2,z3,z4,z5,St,visit.num-1)) 
  # why visit.num-1 not visit.num?
  # reason: visit.num is least 1, but a value of a poisson r.v. can be 0
  countdata2=as.data.frame(cbind(z12,z22,z32,z42,z52,St2,visit.num2-1))
  names(countdata)=c('z1','z2','z3','z4','z5','St','visits')
  names(countdata2)=c('z1','z2','z3','z4','z5','St','visits')
  count <- glm(visits~z1+z2+z3+z4+z5+offset(log(St)),family=poisson(link=log),data=countdata)
  # model 5
  lambdas=exp(as.matrix(cbind(rep(1,ndim),countdata2[,1:5]))%*%as.matrix(count[['coefficients']])[,1])
  # predict the value of Lambda on test dataset (countdata2 here) using the previous model 5 named "count"

N_prect_lis = rep(0, ndim)
for (i in 1:ndim){
  N_prect=ceiling(min(exp_survtime[i],20-X2[i])*lambdas[i])+1
  N_prect_lis[i] = N_prect
}

#model = vglm(SurvS4(visit.num-1, 1 - delta) ~ Z, offset = log(St), cens.poisson)
#summary(model)

#Z2_copy = as.data.frame(cbind(Z2,log(St2)))
#colnames(Z2_copy) = c("z1","z2","z3","z4","z5","log(St)")

#pred_N_train = ceiling(exp(predict(model))) + 1 
#summary(pred_N_train - ceiling(St*lambda))

#ans_1 = cbind(rep(1, dim(Z2)[1]), Z2) %*% as.matrix(model@coefficients)
#pred_CPR = ceiling(St2 * exp(ans_1)) + 1


#dat = data.frame(num = visit.num - 1, ind = 1 - delta, t = St)
#r = (c(as.numeric(model@coefficients),0))  
# X = Z
# neg_log_lik(r, dat, Z)

#fit = optim(f = neg_log_lik, par = r, obj = dat, X = Z, control=list(maxit=10000))
#n_fit = length(fit[["par"]])
#ans = cbind(rep(1, dim(Z2)[1]), Z2) %*% as.matrix(fit[["par"]][-n_fit])
#pred_CGPR = ceiling(St2 * exp(ans)) + 1

true_N = visit.num2
  
  # result = c(mean(N_prect_lis - true_N), mean(pred_CPR - true_N), mean(pred_CGPR - true_N))
  ans = N_prect_lis - true_N
  result = c(min(ans), mean(ans), median(ans), max(ans), sd(ans))
  return(result)
},
error = function(e){return(seed^2)}
)
}
#################################################################################
#  Running the function nSim times:
#################################################################################
nSim = 1000 # test_num
SMNC_PredN_ModelMisSpec <- sfLapply(1:nSim,wrapper)
#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################
save(SMNC_PredN_ModelMisSpec,file = 'SMNC_PredN_ModelMisSpec.RData')
#################################################################################
#  Close all connections:
#################################################################################
sfStop()

