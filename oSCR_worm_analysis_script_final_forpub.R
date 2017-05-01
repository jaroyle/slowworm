load("slowWorms_final.RData")
library(oSCR)
scrFrame <- make.scrFrame(caphist = scrFrame$caphist,
                          traps = scrFrame$traps,
                          indCovs = scrFrame$indCovs,
                          trapCovs = scrFrame$trapCovs,
                          trapOperation = scrFrame$trapOperation)

plot(ssDF,scrFrame,collapse=TRUE)

ss<- rbind(ssDF[[1]],ssDF[[2]],ssDF[[3]])
tiff("fig1.tiff")
  plot(ss,cex=0.5,asp=1)
  points(ssDF[[1]],pch=".",cex=.2,col="gray")
  points(ssDF[[1]],pch=20)
  points(ssDF[[2]],pch=".",cex=.2,col="gray")
  points(ssDF[[2]],pch=20)
  points(ssDF[[3]],pch=".",cex=.2,col="gray")
  points(ssDF[[3]],pch=20)
  text(666600,200100,"area = 1")
  text(666200,199800,"area = 3")
  text(665900, 199450,"area = 2")
dev.off()
 

#sensitivity to state space resolution (3m used in paper)
sesitivity.test <- FALSE
if(sesitivity.test){
  keepers <- matrix(NA,4,3)
  ii <- 1
  for(i in c(40,30,20,10)){
    sss <- make.ssDF(scrFrame,buffer = 70, res = i)
    tmp <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame=scrFrame, ssDF=sss,
                    start.vals=c(-3.6, log(20), -1), multicatch=TRUE, trimS=70)
    keepers[ii,] <- tmp$outStats$mle
    ii <- ii+1
  }
}

#fit the models
new1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame=scrFrame, ssDF=ssDF,
                 start.vals=c(-3.6, log(20), -1), multicatch=TRUE, trimS=70)

new2 <- oSCR.fit(model=list(D~1,p0~session,sig~1), scrFrame=scrFrame, ssDF=ssDF,
                 start.vals=c(-3.6, 0, 0, log(20), -1), multicatch=TRUE, trimS=70)

new3 <- oSCR.fit(model=list(D~session,p0~1,sig~1), scrFrame, ssDF=ssDF, 
                 start.vals=c(-3.6, log(13), -3, 0, 0) , multicatch=TRUE, trimS=70)

new4 <- oSCR.fit(model=list(D~session,p0~session,sig~1), scrFrame, ssDF=ssDF,
                 start.vals=c(-2.9, -0.5, -1, log(13), -3, 0, 0), multicatch=TRUE, trimS=70 )

new5 <- oSCR.fit(model=list(D~session,p0~session + doy,sig~1), scrFrame, ssDF=ssDF,
                 start.vals=c(-2.9, -0.5, -1, log(13),0, -4, 1.2, 1.2), multicatch=TRUE, trimS=70)

new6 <- oSCR.fit(model=list(D~session,p0~session + doy + I(doy^2),sig~1), scrFrame, ssDF=ssDF,
                 start.vals=c(-2.4, -1.2, -1.9, log(13),0,0, -4, 1.2, 1.2), multicatch=TRUE , trimS=70)

new7 <- oSCR.fit(model=list(D~session,p0~1 + doy,sig~1), scrFrame, ssDF=ssDF,
                 start.vals=c(-3.4, log(13),0, -4, 1.2, 1.2), multicatch=TRUE, trimS=70 )

new8 <- oSCR.fit(model=list(D~session,p0~1 + doy + I(doy^2),sig~1), scrFrame, ssDF=ssDF,
                 start.vals=c(-3.4, log(13),0,0, -4, 0.8, 0.5), multicatch=TRUE, trimS=70 )

fl <- fitList.oSCR(list(new1, new2, new3, new4, new5, new6, new7, new8), rename=TRUE)
ms <- modSel.oSCR(fl)

save(new1, new2, new3, new4, new5, new6, new7, new8, file="finalmodelsCS.RData")

library(msm)
mn <- new6$rawOutput$estimate[1:3]
cov <- solve(new6$rawOutput$hessian[1:3, 1:3])
c(plogis(mn[1]), plogis(mn[1]+mn[2]), plogis(mn[1]+mn[3]) )
deltamethod(list(~exp(x1)/(1+exp(x1)), ~exp(x1+x2)/(1+exp(x1+x2)), ~exp(x1+x3)/(1+exp(x1+x3))),
            mean=mn, cov=cov)

mn<-new6$rawOutput$estimate[7:9]
cov<- solve(new6$rawOutput$hessian[7:9, 7:9])
c( exp(mn[1]), exp(mn[1]+mn[2]), exp(mn[1]+mn[3]) )
deltamethod( list(~exp(x1)*10000/9, ~exp(x1+x2)*10000/9, ~exp(x1+x3)*10000/9 ),
 mean=mn, cov=cov)


tiff("fig2.tiff")
X<- cbind(rep(1, length(doy)),doy,doy^2)
parms<- new6$outStats[c(1,5:6),3]
p.fit<- plogis(X%*%parms)
parms<- new6$outStats[c(1,5:6),3]
parms[1]<- parms[1]  + new6$outStats[2,3]
p2.fit<- plogis(X%*%parms)
parms<- new6$outStats[c(1,5:6),3]
parms[1]<- parms[1]  + new6$outStats[3,3]
p3.fit<- plogis(X%*%parms)
xgr<- 7*doy + 60
plot(7*doy + 60,p.fit, type="l",xlab="Day of year", ylab="Detection probability",ylim=c(0,0.10))
lines(xgr, p2.fit, lty=2)
lines(xgr, p3.fit, lty=3)

day.opt<-   (-parms[2]/(2*parms[3]))*7 + 60
lines(c(day.opt,day.opt), c(0,1))
#text(70, 0.005, labels=" optimal")
legend(95,0.10, legend=c("area 1", "area 2", "area 3"), lty=1:3)
dev.off()




pix.area<- 9
# Table in MS:
tab<- ms$coef.tab
tmp<- cbind(0, tab[, 7:8])
tmp[is.na(tmp)]<- 0
D<- exp(tab[,2] + tmp)*(10000/pix.area)
tab<-cbind(D, sigma = exp(tab[,4]), ms$aic.tab[,c("AIC","dAIC")] )
dimnames(tab)<- list(as.character(ms$aic.tab[,1]), c("D1","D2","D3","sigma"
,"AIC","dAIC"))

 sigma<- exp(new6$outStats[4,3] )
 a95<- pi*( sqrt(5.99)*sigma )^2

 
 est<- new6$rawOutput$estimate
 vc<- new6$rawOutput$hessian
 vc<- solve(vc)

# Delta rule approximation for SE of Nhat:
library(msm)
#  SE of area
deltamethod(~   pi*(sqrt(5.99)*exp(x4))^2, est, vc) 




 


