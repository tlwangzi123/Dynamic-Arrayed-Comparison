load("SMNC_MultiNorm.rda")
load("SMNC_Normal_N_CPR.rda")
load("SMNC_Permu_Normal_N.rda")
load("SMNC_Permu_Normal_N_CPR.rda")


v1=r1=v2=r2=NULL
for(i in 1:length(SMNC_MultiNorm)){
  v1=cbind(v1,SMNC_MultiNorm[[i]][,4])
  r1=cbind(r1,SMNC_permu_Normal_N[[i]][,4])
  v2=cbind(v2,SMNC_Normal_N_CPR[[i]][,4])
  r2=cbind(r2,SMNC_permu_Normal_N_CPR[[i]][,4])
}
v1=rowMeans(v1)
r1=rowMeans(r1)
v2=rowMeans(v2)
r2=rowMeans(r2)

#par(mfrow=c(1,1))
#m <- matrix(c(1,2,3,4,5,6),nrow = 2,ncol = 3,byrow = TRUE)
#layout(mat = m,widths = c(0.4,0.4,0.2))
#par(mfrow=c(2,2))

plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,v1,v2,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2)
lines(v2~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=1,col="grey65",lwd=2)
lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col="grey65",lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)


legend("topleft",cex=1.2,legend=c('chisq','chisq using CPR','permu','permu using CPR','Family-wise error = alpha'), col=c("grey0","grey0","grey65","grey65","grey40"), inset = 0, lty=c(1,2,1,2,4),pch=c(NA,NA,15,15,NA),lwd=2)


plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,v1,v2,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)


legend("topleft",cex=1.2,legend=c('chisq','permu','Family-wise error = alpha'), col=c("grey0","grey0","grey40"), inset = 0, lty=c(1,2,4),pch=c(NA,NA,NA),lwd=2)
