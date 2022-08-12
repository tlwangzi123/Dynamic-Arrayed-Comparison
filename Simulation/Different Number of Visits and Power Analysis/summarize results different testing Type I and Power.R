load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_typeI_halfQ.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_typeI_doubleQ.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_typeI_cen30.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_typeI_cen10.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_adjusted.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_power_halfQ.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_power_doubleQ.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_power_cen30.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_power_cen10.rda")
load("C:\\Users\\Zander Wang\\Box\\ZhengWang\\Prospective MNC\\SMNC_permu_MultiT5_power.rda")


r1=r2=r3=r4=r5=NULL
for(i in 1:1000){
  r1=cbind(r1,SMNC_permu_MultiT5_adjusted[[i]][,3])
  r2=cbind(r2,SMNC_permu_MultiT5_typeI_cen10[[i]][,3])
  r3=cbind(r3,SMNC_permu_MultiT5_typeI_cen30[[i]][,3])
  r4=cbind(r4,SMNC_permu_MultiT5_typeI_halfQ[[i]][,3])
  r5=cbind(r5,SMNC_permu_MultiT5_typeI_doubleQ[[i]][,3])
}
r1=rowMeans(r1)
r2=rowMeans(r2)
r3=rowMeans(r3)
r4=rowMeans(r4)
r5=rowMeans(r5)

par(mfrow=c(1,2))
plot(r1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,r4,r5)), col='gray0' ,cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Under Null')

lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col='gray20',lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col='gray40',lwd=2)
lines(r4~seq(0.001,0.1,by=0.0001),lty=4,col='gray60',lwd=2)
lines(r5~seq(0.001,0.1,by=0.0001),lty=5,col='gray80',lwd=2)
legend("bottomright",legend=c("(1) Original","(2) Shorter Study","(3) Longer Study",
                              "(4) Less Visits","(5) More Visits"),
       col=c(col='gray0',col='gray20',col='gray40',col='gray60',col='gray80'),lwd=3,lty=1:5,pch=rep(NA,5))



r1=r2=r3=r4=r5=NULL
for(i in 1:1000){
  r1=cbind(r1,SMNC_permu_MultiT5_power[[i]])
  r2=cbind(r2,SMNC_permu_MultiT5_power_cen10[[i]])
  r3=cbind(r3,SMNC_permu_MultiT5_power_cen30[[i]])
  r4=cbind(r4,SMNC_permu_MultiT5_power_halfQ[[i]])
  r5=cbind(r5,SMNC_permu_MultiT5_power_doubleQ[[i]])
}
r1=rowMeans(r1)
r2=rowMeans(r2)
r3=rowMeans(r3)
r4=rowMeans(r4)
r5=rowMeans(r5)

plot(r1~seq(0.001,0.1,by=0.0001),ylim=c(0,max(r1,r2,r3,r4,r5)), col='gray0' ,cex.main=1.5, cex.lab=1.3,type='l',xlab='alpha',ylab='Family-wise Error Rate',lwd=2,main='Under Alternative')

lines(r2~seq(0.001,0.1,by=0.0001),lty=2,col='gray20',lwd=2)
lines(r3~seq(0.001,0.1,by=0.0001),lty=3,col='gray40',lwd=2)
lines(r4~seq(0.001,0.1,by=0.0001),lty=4,col='gray60',lwd=2)
lines(r5~seq(0.001,0.1,by=0.0001),lty=5,col='gray80',lwd=2)
legend("bottomright",legend=c("(1) Original","(2) Shorter Study","(3) Longer Study",
                              "(4) Less Visits","(5) More Visits"),
       col=c(col='gray0',col='gray20',col='gray40',col='gray60',col='gray80'),lwd=3,lty=1:5,pch=rep(NA,5))













