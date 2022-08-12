load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_Gamma4.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_Gamma25.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_Gamma100.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_Gamma4.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_Gamma25.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_Gamma100.rda")

v1=r1=v2=r2=v3=r3=NULL
for(i in 1:length(SMNC_Gamma4)){
  v1=cbind(v1,SMNC_Gamma4[[i]][,1])
  r1=cbind(r1,SMNC_permu_Gamma4[[i]][,1])
  v2=cbind(v2,SMNC_Gamma25[[i]][,1])
  r2=cbind(r2,SMNC_permu_Gamma25[[i]][,1])
  v3=cbind(v3,SMNC_Gamma100[[i]][,1])
  r3=cbind(r3,SMNC_permu_Gamma100[[i]][,1])
}
v1=rowMeans(v1)
r1=rowMeans(r1)
v2=rowMeans(v2)
r2=rowMeans(r2)
v3=rowMeans(v3)
r3=rowMeans(r3)

m <- matrix(c(1,2,3,4,5,6),nrow = 2,ncol = 3,byrow = TRUE)

layout(mat = m,widths = c(0.4,0.4,0.2))

plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,v1,v2,v3,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Constant Mean Trend')
lines(v2~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(v3~seq(0.001,0.1,by=0.0001),lty=3,lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=1,col="grey65",lwd=2)
lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col="grey65",lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col="grey65",lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)

v1=r1=v2=r2=v3=r3=NULL
for(i in 1:length(SMNC_Gamma4)){
  v1=cbind(v1,SMNC_Gamma4[[i]][,2])
  r1=cbind(r1,SMNC_permu_Gamma4[[i]][,2])
  v2=cbind(v2,SMNC_Gamma25[[i]][,2])
  r2=cbind(r2,SMNC_permu_Gamma25[[i]][,2])
  v3=cbind(v3,SMNC_Gamma100[[i]][,2])
  r3=cbind(r3,SMNC_permu_Gamma100[[i]][,2])
}
v1=rowMeans(v1)
r1=rowMeans(r1)
v2=rowMeans(v2)
r2=rowMeans(r2)
v3=rowMeans(v3)
r3=rowMeans(r3)

plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,v1,v2,v3,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Linear Mean Trend')
lines(v2~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(v3~seq(0.001,0.1,by=0.0001),lty=3,lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=1,col="grey65",lwd=2)
lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col="grey65",lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col="grey65",lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",cex=1.2,legend=c('shape 4, scale 5','shape 25, scale 2','shape 100, scale 1','chisq test','permutation test','Family-wise error = alpha'),
       col=c("grey0","grey0","grey0","grey0","grey65","grey40"), inset = 0, lty=c(1,2,3,NA,NA,4),pch=c(NA,NA,NA,15,15,NA),lwd=2)


par(mar= c(5, 4, 4, 2) + 0.1)

v1=r1=v2=r2=v3=r3=NULL
for(i in 1:length(SMNC_Gamma4)){
  v1=cbind(v1,SMNC_Gamma4[[i]][,3])
  r1=cbind(r1,SMNC_permu_Gamma4[[i]][,3])
  v2=cbind(v2,SMNC_Gamma25[[i]][,3])
  r2=cbind(r2,SMNC_permu_Gamma25[[i]][,3])
  v3=cbind(v3,SMNC_Gamma100[[i]][,3])
  r3=cbind(r3,SMNC_permu_Gamma100[[i]][,3])
}
v1=rowMeans(v1)
r1=rowMeans(r1)
v2=rowMeans(v2)
r2=rowMeans(r2)
v3=rowMeans(v3)
r3=rowMeans(r3)

plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,v1,v2,v3,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Quadratic Mean Trend')
lines(v2~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(v3~seq(0.001,0.1,by=0.0001),lty=3,lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=1,col="grey65",lwd=2)
lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col="grey65",lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col="grey65",lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)

v1=r1=v2=r2=v3=r3=NULL
for(i in 1:length(SMNC_Gamma4)){
  v1=cbind(v1,SMNC_Gamma4[[i]][,4])
  r1=cbind(r1,SMNC_permu_Gamma4[[i]][,4])
  v2=cbind(v2,SMNC_Gamma25[[i]][,4])
  r2=cbind(r2,SMNC_permu_Gamma25[[i]][,4])
  v3=cbind(v3,SMNC_Gamma100[[i]][,4])
  r3=cbind(r3,SMNC_permu_Gamma100[[i]][,4])
}
v1=rowMeans(v1)
r1=rowMeans(r1)
v2=rowMeans(v2)
r2=rowMeans(r2)
v3=rowMeans(v3)
r3=rowMeans(r3)

plot(v1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,v1,v2,v3,0.1)),cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Cubic Mean Trend')
lines(v2~seq(0.001,0.1,by=0.0001),lty=2,lwd=2)
lines(v3~seq(0.001,0.1,by=0.0001),lty=3,lwd=2)
lines(r1~seq(0.001,0.1,by=0.0001),lty=1,col="grey65",lwd=2)
lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col="grey65",lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col="grey65",lwd=2)
lines(c(0.001,.1),c(0.001,.1),lty=4,col="grey40",lwd=2)

