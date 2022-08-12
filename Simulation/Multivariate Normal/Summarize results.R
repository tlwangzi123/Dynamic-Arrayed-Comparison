load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_MultiNorm.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_MultiNorm2.rda")
par(mfrow=c(2,2))
r1=r2=r3=r4=NULL
for(i in 1:length(SMNC_MultiNorm)){
  r1=cbind(r1,SMNC_MultiNorm[[i]][,1])
  r2=cbind(r2,SMNC_MultiNorm[[i]][,1])
  r3=cbind(r3,SMNC_MultiNorm[[i]][,1])
  r4=cbind(r4,SMNC_MultiNorm[[i]][,1])
}
r1=rowMeans(r1)
r2=rowMeans(r2)
r3=rowMeans(r3)
r4=rowMeans(r4)
par(mfrow=c(2,2))
plot(r1~seq(0.001,0.1,by=0.0001),cex.main=1.5, cex.lab=1.3,ylim=c(0,0.1),type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Constant Mean Trend')
lines(c(0.001,.1),c(0.001,.1),lty=2,col="grey50",lwd=2)

plot(r2~seq(0.001,0.1,by=0.0001),cex.main=1.5, cex.lab=1.3,ylim=c(0,0.1),type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Linear Mean Trend')
lines(c(0.001,.1),c(0.001,.1),lty=2,col="grey50",lwd=2)

plot(r3~seq(0.001,0.1,by=0.0001),cex.main=1.5, cex.lab=1.3,ylim=c(0,0.1),type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Quadratic Mean Trend')
lines(c(0.001,.1),c(0.001,.1),lty=2,col="grey50",lwd=2)

plot(r4~seq(0.001,0.1,by=0.0001),cex.main=1.5, cex.lab=1.3,ylim=c(0,0.1),type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Cubic Mean Trend')
lines(c(0.001,.1),c(0.001,.1),lty=2,col="grey50",lwd=2)
legend("bottomright",legend=c('From simulation','Family-wise error = alpha'),
       col=c("grey0","grey50"),lty=c(1,2),pch=c(NA,NA),lwd=2)

