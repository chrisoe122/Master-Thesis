source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Final code/SGBM-fun.R")
source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Final code/LSM-fun.R")
source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Final code/CHE-fun-Glau.R")

library(scales)



### Generating data

#1x5Y

#Parameters 
fm<-0.01
KQ<-100*10^3
KA<-100*10^3
TT<-100

dt<-0.05

X5<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
Xa5<- matrix(data=NA, nrow=KA, ncol=(TT+1))


set.seed(2023)

#Under Q
sigma<-0.02
lambda<-0.02

#Generating X under Q
X5[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  X5[,i]<- r((i-1)*dt,i*dt,X5[,i-1],KQ)
}

#Under P
sigma<-0.01 
lambda<-0.015 

#Generating X under P
Xa5[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  Xa5[,i]<- r((i-1)*dt,i*dt,Xa5[,i-1],KA)
}

#Algorithm
CHE51<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4, kk=5, N=80, exercisestep = 20)
CHE52<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1, kk=5, N=80, exercisestep = 20)
CHE53<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6, kk=5, N=80, exercisestep = 20)



SGBM51<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4)
SGBM52<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1)
SGBM53<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6)

LSM51<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4)
LSM52<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1)
LSM53<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6)


###  4x10Y  #####

#Parameters 
fm<-0.01
KQ<-100*10^3
KA<-100*10^3
TT<-200
dt<-0.05

X10<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
Xa10<- matrix(data=NA, nrow=KA, ncol=(TT+1))

set.seed(2023)


#Under Q
sigma<-0.01
lambda<-0.012

#Generating X under Q 
X10[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  X10[,i]<- r((i-1)*dt,i*dt,X10[,i-1],KQ)
}

#Under P
sigma<-0.006 
lambda<-0.008

#Generating X under P
Xa10[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  Xa10[,i]<- r((i-1)*dt,i*dt,Xa10[,i-1],KA)
}

#Algorithm
CHE101<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4, kk=5, N=80, exercisestep=20)
CHE102<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep=20) 
CHE103<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1.6, kk=5, N=80, exercisestep=20)



SGBM101<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4)
SGBM102<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
SGBM103<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1.6)

LSM101<-LSM(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4)
LSM102<-LSM(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)  
LSM103<-LSM(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1.6)


######## TABLE VALUES ########

#1x5Y
CHE51[1] 
CHE52[1]
CHE53[1]

SGBM51[1] 
SGBM52[1]
SGBM53[1]

LSM51[1]
LSM52[1]
LSM53[1]

#4x10Y
CHE101[1] 
CHE102[1]
CHE103[1]

SGBM101[1]
SGBM102[1]
SGBM103[1]

LSM101[1]
LSM102[1]
LSM103[1]



######  Tid  ##############
#1x5Y
CHE51[9] 
CHE52[9]
CHE53[9]

SGBM51[9] 
SGBM52[9]
SGBM53[9]



LSM51[9]
LSM52[9]
LSM53[9]

#4x10Y
CHE101[9] 
CHE102[9]
CHE103[9]

SGBM101[9]
SGBM102[9]
SGBM103[9]

LSM101[9]
LSM102[9]
LSM103[9]


########TJEK AF UDE AF DOMAIN
CHE102[11]
PO<-sort(X10, decreasing = T)
PO[1:10]
max(X10[,190])
#1
#3
#5
#6
#4
#2
#7
#9
#8

#Done after 10. Then all in the domain forever
EECHE2<-CHE102[[3]]
EECHE1<-CHE101[[3]]
EECHE3<-CHE103[[3]]
Xsort<-X10[order(X10[,201], decreasing = T),]
EECHE1<-EECHE1[order(X10[,201], decreasing = T),]
EECHE2<-EECHE2[order(X10[,201], decreasing = T),]
EECHE3<-EECHE3[order(X10[,201], decreasing = T),]
Xsort[1,]
EECHE1[1,]
EECHE2[1,]
EECHE3[1,]
#Exercise in all algorithms.






######### PLOTS ##############
names(SGBM102)
windows(width=10, height=8)


plot(seq(0,10,0.05),SGBM102[[7]], type="l", ylim=c(0,7), lwd=2, main="", ylab="EE", xlab = "time")
grid(nx=NULL)
lines(seq(0,10,0.05),SGBM102[[8]], type="l", lty=2, lwd=2)
lines(seq(0,10,0.05),LSM102[[7]], type="l", lty=1, lwd=2, col="blue")
lines(seq(0,10,0.05),LSM102[[8]], type="l", lty=2, lwd=2,col="blue")
lines(seq(0,10,0.05),CHE102[[7]], type="l", lty=1, lwd=2, col="red")
lines(seq(0,10,0.05),CHE102[[8]], type="l", lty=2, lwd=2, col="red")
legend(x="bottomleft", legend=c("SGBM", "LSM", "CHE"), lwd=2, col=c("black", "blue", "red"))



###PFE PLOT
plot(seq(0,10,0.05),SGBM102[[5]], type="l", ylim=c(0,35), lwd=2, main="", ylab="PFE", xlab = "time")
grid(nx=NULL)
lines(seq(0,10,0.05),SGBM102[[6]], type="l", lty=2, lwd=2)
lines(seq(0,10,0.05),LSM102[[6]], type="l", lty=2, lwd=2, col="blue")
lines(seq(0,10,0.05),LSM102[[5]], type="l", lwd=2, col="blue")
lines(seq(0,10,0.05),CHE102[[6]], type="l", lty=2, lwd=2, col="red")
lines(seq(0,10,0.05),CHE102[[5]], type="l", lwd=2, col="red")
legend(x="bottomleft", legend=c("SGBM", "LSM", "CHE"), lwd=2, col=c("black", "blue", "red"))


#### Fig 4a
plot(seq(2,3.95,0.05),SGBM102[[5]][41:80], type="l", ylim=c(14,20), lwd=2, main="", ylab="PFE", xlab = "time")
grid(nx=NULL)
lines(seq(2,3.95,0.05),LSM102[[5]][41:80], type="l", lwd=2, col="blue")
lines(seq(2,3.95,0.05),CHE102[[5]][41:80], type="l", lwd=2, col="red", lty=2)
legend(x="topleft", legend=c("SGBM", "LSM", "CHE"), lwd=2, lty=c(1,1,2),  col=c("black", "red", "blue"))


#### Fig 4b

plot(seq(6,7.95,0.05),SGBM102[[5]][121:160], ylim=c(2,5.5), type="l", lwd=2, main="", ylab="PFE", xlab = "time")
grid(nx=NULL)
lines(seq(6,7.95,0.05),LSM102[[5]][121:160], type="l", lwd=2, col="blue")
lines(seq(6,7.95,0.05),CHE102[[5]][121:160], type="l", lwd=2, lty=2,col="red")
legend(x="bottomleft", legend=c("SGBM", "LSM", "CHE"), lwd=2, lty=c(1,1,2), col=c("black", "blue", "red"))


### Fig 5a
step<-126
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, ylim=c(0,70),main="", ylab="Continuation value", xlab = "Interest rate")
CHEord<-CHE102[[2]][order(X10[,step]),]
lines(Xord[,step],CHEord[,step], col="red", lty=2,lwd=2)
CLSM<-LSM102[[10]][1,step]+LSM102[[10]][2,step]*LSM102[[11]][,step]+LSM102[[10]][3,step]*LSM102[[12]][,step]+
  LSM102[[10]][4,step]*LSM102[[13]][,step]
CLSMord<- CLSM[order(X10[,step])]
lines(Xord[,step],CLSMord, col="blue", lwd=2)
grid(nx=NULL)
legend(x="bottomleft", legend=c("SGBM", "LSM", "CHE"), lwd=2, lty=c(1,1,2),col=c("black", "blue", "red"))

### Fig 5b
qd<-density(X10[,126])
ad<-density(Xa10[,126])

plot(qd, ylim=c(0,30), type='l', lwd=2, main="", ylab="density", xlab = "Interest rate")
grid(nx=NULL)
lines(ad, col='red', lwd=2)
legend(x="topleft", legend=c("q-density", "p-density"), lwd=2, col=c("black", "red"))



#################################### ANDRE PLOTS ##################################################

### BEREGNING AF CUT-OFF POINT. TIME=9
CutX<-min(X10[(SGBM102[[10]][,181]==SGBM102[[2]][,181]),181])
CutY<-SGBM102[[10]][X10[,181]==CutX,181]
mean(CutY>=SGBM102[[2]][(SGBM102[[10]][,181]==SGBM102[[2]][,181]),181])# =1
mean(CutY<SGBM102[[2]][(SGBM102[[10]][,181]==SGBM102[[2]][,181]),181])# = 0 #Så det er legit cut-off point

CHE20<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=20, exercisestep=20)
############ VALUE OF THE OPTION AT EXERCISE TIME AND BEFORE. ZOOM AND WHOLE

windows(width=10, height=8)

step<-181
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, main="", ylab="Value", xlab = "Interest rate")
#VLSM<-LSM102[[2]][order(X10[,step]),]
#points(Xord[,step],VLSM[,step], col=rgb(red = 1, green = 0, blue = 0, alpha = 0.1), cex=1)
CHEord<-CHE102[[2]][order(X10[,step]),]
CHEord1<-CHE20[[2]][order(X10[,step]),]
CHEord2<-CHE40[[2]][order(X10[,step]),]
CHEord3<-CHE110[[2]][order(X10[,step]),]
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
lines(Xord[,step],CHEord1[,step], col="blue", lwd=2, lty=2)
lines(Xord[,step],CHEord2[,step], col="green", lwd=2, lty=2)
lines(Xord[,step],CHEord3[,step], col="purple", lwd=2, lty=2)
grid(nx=NULL)
legend(x="bottomleft", legend=c("SGBM", "CHE20", "CHE40","CHE80","CHE110"), col=c("black", "blue","green", "orange", "purple"), lty=c(1,2,2,2,2), lwd=2, bg="white")

step<-180
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, main="", ylab="Value", xlab = "Interest rate")
#VLSM<-LSM102[[2]][order(X10[,step]),]
#points(Xord[,step],VLSM[,step], col=rgb(red = 1, green = 0, blue = 0, alpha = 0.1), cex=1)
CHEord<-CHE102[[2]][order(X10[,step]),]
CHEord1<-CHE20[[2]][order(X10[,step]),]
CHEord2<-CHE40[[2]][order(X10[,step]),]
CHEord3<-CHE110[[2]][order(X10[,step]),]
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
lines(Xord[,step],CHEord1[,step], col="blue", lwd=2, lty=2)
lines(Xord[,step],CHEord2[,step], col="green", lwd=2, lty=2)
lines(Xord[,step],CHEord3[,step], col="purple", lwd=2, lty=2)
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
grid(nx=NULL)
legend(x="bottomleft", legend=c("SGBM", "CHE20", "CHE40","CHE80","CHE110"), col=c("black", "blue","green", "orange", "purple"), lty=c(1,2,2,2,2), lwd=2, bg="white")





################# EKSTRA: LAV PÆÆææææææææææææææææææNT. ##############
plot(X10[,step],LSM102[[2]][,step])
Great<- (X10[,step]>0.03)
Und<- (X10[,step]<0.03005)
Her<-as.logical(Great*Und)
LSM102[[2]][Her,step]
X10[Her,201]
U(0.0092,5,5,0.0109,0.01,0.012)
mean(LSM102[[2]][,step]==1.79660771)
mean(Her)
plot(LSM102[[2]][100,])




step<-180
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, main="Value before exercise time (t=8.95)", ylab="Value", xlab = "Interest rate")
CHEord<-CHE102[[2]][order(X10[,step]),]
grid(nx=NULL)

CLSMd<-matrix(data=NA, nrow=KQ, ncol=4)
CLSMd[,1]<-1
CLSMd[,2]<-X10[,step]

for (i in 3:4){
  CLSMd[,i]<-CLSMd[,i-1]*CLSMd[,2]
}
CLSM<-CLSMd%*%LSM102[[10]][,step]
CLSMord<- CLSM[order(X10[,step])]
points(Xord[,step],CLSMord, col="red", lwd=2)
lines(Xord[,step],Cmord[,step], lwd=2)
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
legend(x="bottomleft", legend=c("SGBM", "CHE", "LSM"), col=c("black", "orange", "red"), lty=c(1,1,0), lwd=2, pch=c(NA,NA,1), bg="white")


#ZOOM
step<-181
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, ylim=c(0,1),xlim=c(0.0,0.05), main="", ylab="Value", xlab = "Interest rate")
CHEord<-CHE102[[2]][order(X10[,step]),]
CHEord1<-CHE20[[2]][order(X10[,step]),]
CHEord2<-CHE40[[2]][order(X10[,step]),]
CHEord3<-CHE110[[2]][order(X10[,step]),]
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
lines(Xord[,step],CHEord1[,step], col="blue", lwd=2, lty=2)
lines(Xord[,step],CHEord2[,step], col="green", lwd=2, lty=2)
lines(Xord[,step],CHEord3[,step], col="purple", lwd=2, lty=2)


grid(nx=NULL)
points(CutX,CutY, col="red", pch=16, cex=1.1)
legend(x="bottomleft", legend=c("SGBM", "CHE20", "CHE40","CHE80","CHE110"), col=c("black", "orange","blue", "green", "purple"), lty=c(1,2,2,2,2), lwd=2, bg="white")




step<-181
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]
plot(Xord[,step],Cmord[,step], type="l",lwd=2, ylim=c(0,1),xlim=c(0.0,0.05), main="", ylab="Value", xlab = "Interest rate")
CHEord<-CHE102[[2]][order(X10[,step]),]
CHEord1<-CHE20[[2]][order(X10[,step]),]
CHEord2<-CHE40[[2]][order(X10[,step]),]
CHEord3<-CHE110[[2]][order(X10[,step]),]
lines(Xord[,step],CHEord[,step], col="orange", lwd=2, lty=2)
lines(Xord[,step],CHEord1[,step], col="blue", lwd=2, lty=2)
lines(Xord[,step],CHEord2[,step], col="green", lwd=2, lty=2)
lines(Xord[,step],CHEord3[,step], col="purple", lwd=2, lty=2)
grid(nx=NULL)
legend(x="bottomleft", legend=c("SGBM", "CHE20", "CHE40","CHE80","CHE110"), col=c("black", "blue","green", "orange", "purple"), lty=c(1,2,2,2,2), lwd=2, bg="white")














############################ SPÆNDENDE PLOTS. SER UD TIL AT DE BLOT ER UENIGE VED EXERCISE DATES. ELLERS FØLGES DE AD.
windows(10,8)

VmCHE<-rep(NA,201)
VmCHE40<-rep(NA,201)
VmCHE110<-rep(NA,201)
VmSGBM<-rep(NA,201)

CHE110<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=20)
CHE40<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=20)


for (i in 1:191){
  VmCHE[i]<-mean(CHE102[[2]][,i])
  VmCHE40[i]<-mean(CHE40[[2]][,i])
  VmCHE110[i]<-mean(CHE110[[2]][,i])
  VmSGBM[i]<-mean(SGBM102[[2]][,i])
}

ab<-182
bc<-201
plot(seq(0,10,0.05)[ab:bc],VmSGBM[ab:bc], ylim=c(1.02,1.07),type="l", lwd=2, main="V from SGBM and CHE", ylab="V", xlab = "time")
grid(nx=NULL)
lines(seq(0,10,0.05)[ab:bc],VmCHE[ab:bc], type="l", lwd=2, col="blue")
#lines(seq(0,10,0.05)[ab:bc],VmCHEvG[ab:bc], type="l", lwd=2, col="red")
legend(x="bottomleft", legend=c("SGBM", "CHE"), col=c("black", "blue"), lty=c(1,1), lwd=3, bg="white")


ab<-162
bc<-181
plot(seq(0,10,0.05)[ab:bc],VmSGBM[ab:bc], type="l", lwd=2, main="", ylab="Value", xlab = "Time")
grid(nx=NULL)
lines(seq(0,10,0.05)[ab:bc],VmCHE[ab:bc], type="l", lwd=2, col="blue")
#lines(seq(0,10,0.05)[ab:bc],VmCHEvG[ab:bc], type="l", lwd=2, col="red")
legend(x="bottomleft", legend=c("SGBM", "CHE"), col=c("black", "blue"), lty=c(1,1), lwd=3, bg="white")


ab<-142
bc<-161
plot(seq(0,10,0.05)[ab:bc],VmSGBM[ab:bc], type="l", lwd=2, main="", ylab="Value", xlab = "Time")
grid(nx=NULL)
lines(seq(0,10,0.05)[ab:bc],VmCHE[ab:bc], type="l", lwd=2, col="blue")
#lines(seq(0,10,0.05)[ab:bc],VmCHEvG[ab:bc], type="l", lwd=2, col="red")
legend(x="bottomleft", legend=c("SGBM", "CHE"), col=c("black", "blue"), lty=c(1,1), lwd=3, bg="white")


ab<-1
bc<-81
plot(seq(0,10,0.05)[ab:bc],VmSGBM[ab:bc], ylim=c(6.12,6.26), type="l", lwd=2, main="", ylab="Value", xlab = "Time")
grid(nx=NULL)
lines(seq(0,10,0.05)[ab:bc],VmCHE[ab:bc], type="l", lwd=2, col="blue")
lines(seq(0,10,0.05)[ab:bc],VmCHE40[ab:bc], type="l", lwd=2, col="red")
lines(seq(0,10,0.05)[ab:bc],VmCHE110[ab:bc], type="l", lwd=2, col="green")
legend(x="topleft", legend=c("SGBM", "CHE110", "CHE80", "CHE40"), col=c("black", "green","blue", "red"), lty=c(1,1), lwd=3, bg="white")





########## FIGUR:"Different between mean of Cont. value of SGBM and CHE (1:195)"

holdd<-rep(NA,191) #De sidste pkt. er lidt ustabile, da max(X) er udenfor interval
for (i in 1:191){
  holdd[i]<-VmSGBM[i]-VmCHE[i]
}

hold40<-rep(NA,191) #De sidste pkt. er meget ustablie.
for(i in 1:191){
  hold40[i]<-VmSGBM[i]-VmCHE40[i]
}

hold110<-rep(NA,191) #De sidste pkt. er meget ustablie.
for(i in 1:191){
  hold110[i]<-VmSGBM[i]-VmCHE110[i]
}

plot(seq(0,9.5,0.05),holdd, col="blue", main="", xlab="Time",ylab="", ylim=c(-0.02,0.08))
points(seq(0,9.5,0.05),hold40, col="red")
points(seq(0,9.5,0.05),hold110)
segments(-5,0,205,0, lty=2)
legend(x="topright", legend=c("CHE110", "CHE80", "CHE40"), col=c("black","blue", "red"), lty=c(1,1), lwd=3, bg="white")

#Plot'sne viser, at der hvor de er mest uenige med hinanden er ved exercising.









############# FIGUR ############
windows(10,8)
plot(X10[,185],CHE102[[2]][,185], main="", ylab="Value", xlab="Interest rate", ylim=c(0,20), xlim=c(-0.16,0.16))
points(CHE102[[11]], CHE102[[10]][,185], col="red", pch=16)
plot(X10[,85],CHE102[[2]][,85], main="", ylab="Value", xlab="Interest rate", ylim=c(0,150), xlim=c(-0.16,0.16))
points(CHE102[[11]], CHE102[[10]][,85], col="red", pch=16)
plot(CHE102[[11]], CHE102[[10]][,199], main="Option valued at Chebyshev points. Time=9.90", ylab="V", xlab="Interest rate")
plot(CHE102[[11]], CHE102[[10]][,121], main="Option valued at Chebyshev points. Time=6", ylab="V", xlab="Interest rate")
plot(CHE102[[11]], CHE102[[10]][,45], main="Option valued at Chebyshev points. Time=2.20", ylab="V", xlab="Interest rate")









######### Creating SD ###################
## 1x4Y
LSMsd51<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd51<-matrix(data=NA, nrow=4, ncol=10)
CHEsd51<-matrix(data=NA, nrow=4, ncol=10)
LSMsd52<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd52<-matrix(data=NA, nrow=4, ncol=10)
CHEsd52<-matrix(data=NA, nrow=4, ncol=10)
LSMsd53<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd53<-matrix(data=NA, nrow=4, ncol=10)
CHEsd53<-matrix(data=NA, nrow=4, ncol=10)


set.seed(2023)
for (w in 1:10){
  #Parameters 
  fm<-0.01
  KQ<-100*10^2
  KA<-100*10^2
  TT<-100
  dt<-0.05
  
  X5<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
  Xa5<- matrix(data=NA, nrow=KA, ncol=(TT+1))
  
  #Under Q
  sigma<-0.02
  lambda<-0.02
  
  #Generating X under Q
  X5[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    X5[,i]<- r((i-1)*dt,i*dt,X5[,i-1],KQ)
  }
  
  #Under P
  sigma<-0.01 
  lambda<-0.015 
  
  #Generating X under P
  Xa5[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    Xa5[,i]<- r((i-1)*dt,i*dt,Xa5[,i-1],KA)
  }
  print(w)
  CHEsd51[,w]<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4, kk=5, N=80, exercisestep = 20)[[1]]
  CHEsd52[,w]<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1, kk=5, N=80, exercisestep = 20)[[1]]
  CHEsd53[,w]<-CheVGLAU(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6, kk=5, N=80, exercisestep = 20)[[1]]
  
  SGBMsd51[,w]<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4)[[1]]
  SGBMsd52[,w]<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1)[[1]]
  SGBMsd53[,w]<-SGBM(X=X5, Xa=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6)[[1]]
  
  LSMsd51[,w]<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=0.4)[[1]]
  LSMsd52[,w]<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1)[[1]]
  LSMsd53[,w]<-LSM(X=X5, XaL=Xa5, KQ=KQ, KA=KA, dt=dt, Q=1, EX=4, fm=fm, sigma=0.02, lambda=0.02, kp=1.6)[[1]]
  
  print(w)
}

apply(CHEsd51,1,sd)
apply(CHEsd52,1,sd)
apply(CHEsd53,1,sd)

apply(SGBMsd51,1,sd)
apply(SGBMsd52,1,sd)
apply(SGBMsd53,1,sd)

apply(LSMsd51,1,sd)
apply(LSMsd52,1,sd)
apply(LSMsd53,1,sd)




#4x10Y
LSMsd101<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd101<-matrix(data=NA, nrow=4, ncol=10)
CHEsd101<-matrix(data=NA, nrow=4, ncol=10)
LSMsd102<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd102<-matrix(data=NA, nrow=4, ncol=10)
CHEsd102<-matrix(data=NA, nrow=4, ncol=10)
LSMsd103<-matrix(data=NA, nrow=4, ncol=10)
SGBMsd103<-matrix(data=NA, nrow=4, ncol=10)
CHEsd103<-matrix(data=NA, nrow=4, ncol=10)


set.seed(2023)
for (w in 1:10){
  #Parameters 
  fm<-0.01
  KQ<-100*10^3
  KA<-100*10^3
  TT<-200
  dt<-0.05
  
  X10<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
  Xa10<- matrix(data=NA, nrow=KA, ncol=(TT+1))
  
  #Under Q
  sigma<-0.01
  lambda<-0.012
  
  #Generating X under Q
  X10[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    X10[,i]<- r((i-1)*dt,i*dt,X10[,i-1],KQ)
  }
  
  #Under P
  sigma<-0.006
  lambda<-0.008 
  
  #Generating X under P
  Xa10[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    Xa10[,i]<- r((i-1)*dt,i*dt,Xa10[,i-1],KA)
  }
  print(w)
  #CHEsd101[,w]<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4, kk=5, N=80, exercisestep = 20)[[1]]
  #CHEsd102[,w]<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep = 20)[[1]]
  #CHEsd103[,w]<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1.6, kk=5, N=80, exercisestep = 20)[[1]]
  
  #SGBMsd101[,w]<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4)[[1]]
  #SGBMsd102[,w]<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]]
  #SGBMsd103[,w]<-SGBM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1.6)[[1]]
  
  #LSMsd101[,w]<-LSM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=0.4)[[1]]
  LSMsd102[,w]<-LSM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]]
  #LSMsd103[,w]<-LSM(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]]
  
  print(w)
}

apply(CHEsd101,1,sd)
apply(CHEsd102,1,sd)
apply(CHEsd103,1,sd)

apply(SGBMsd101,1,sd)
apply(SGBMsd102,1,sd)
apply(SGBMsd103,1,sd)

apply(LSMsd101,1,sd)
apply(LSMsd102,1,sd)
apply(LSMsd103,1,sd)


mean(LSMsd102[1,])
mean(LSMsd103[1,])

LSMsd102[1,]
LSMsd103[4,]

stripchart(LSMsd102[1,])
stripchart(LSMsd102[1,], vertical=T, ylim=c(6.15,6.25), pch=1, col="red")
abline(6.199,0)




############################## Comparision when dt gets smaller
TT<-100
dt<-0.1

X101<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
Xa101<- matrix(data=NA, nrow=KA, ncol=(TT+1))

set.seed(2023)


#Under Q
sigma<-0.01
lambda<-0.012

#Generating X under Q 
X101[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  X101[,i]<- r((i-1)*dt,i*dt,X101[,i-1],KQ)
}

#Under P
sigma<-0.006 
lambda<-0.008

#Generating X under P
Xa101[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  Xa101[,i]<- r((i-1)*dt,i*dt,Xa101[,i-1],KA)
}

SGBM102NYTT<-SGBM(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
SGBM102NYTT<-SGBM(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
SGBM102NYT<-SGBM(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)


###LSM
#poly(3)
LSM102NYTT<-LSM(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM102NYTT<-LSM(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM102NYTT[9]
LSM102NYT<-LSM(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)

#poly(8)
LSM8_10<-LSM8(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM8_20<-LSM8(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM8_40<-LSM8(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM8_100<-LSM8(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)


#poly(5)
LSM5_10<-LSM15(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM5_20<-LSM15(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM5_40<-LSM15(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM5_100<-LSM15(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)


###CHE
# n=110
CHE102NYT<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=10)
CHE102NYT1<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=20)
CHE102NYT<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=50)
CHE102NYTT<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=100)


# n=80
CHE80_10<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep=10)
CHE80_20<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep=20)
CHE80_50<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep=50)
CHE80_100<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=80, exercisestep=100)

# n=40
CHE40_10<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.1, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=10)
CHE40_20<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=20)
CHE40_50<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.02, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=50)
CHE40_100<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.01, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=100)


SGBM102NYT[9]
SGBM102NYTT[9]
CHE102NYT[1]
CHE102NYT[9]
LSM102NYT[1]
LSM102NYTT[1]

CHE80_20[9]
CHE40_10[1]
CHE80_50[9]
CHE40_50[9]
CHE80_100[9]
CHE40_100[9]

CHE80_10[9]
CHE40_20[1]
CHE80_50[1]
CHE40_50[1]
CHE80_100[9]
CHE40_100[1]

LSM8_10[1]
LSM5_20[1]
LSM8_40[9]
LSM5_50[1]
LSM8_100[1]
LSM5_100[1]

LSM5_10[9]
LSM5_40[1]
LSM8_50[9]
LSM5_50[9]
LSM8_100[9]
LSM5_100[9]



##Plots
windows(10,8)
#CHE V0
CHEval40<-c(6.128,6.128,6.137,6.143)
CHEval80<-c(6.182, 6.185,6.188,6.189)
CHEval128<-c(6.189,6.191,6.195,6.196)
SGBMval<-c(6.199,6.200,6.199,6.199)
xdt<-c(0.1,0.05,0.02,0.01)
plot(xdt, CHEval40, type="b", ylim=c(6.1,6.22), xlab="dt", ylab="V0")
lines(xdt, CHEval80, type="b", col="red")
lines(xdt, CHEval128, type="b", col="blue")
lines(xdt, SGBMval, type="b", col="orange")
legend(x="bottomleft", legend=c("CHE40","CHE80","CHE110","SGBM"), col=c("black", "red", "blue", "orange"), lwd=1,bg="white")


#CHE EPE
CHEval40<-c(2.544,2.569,2.584,2.584)
CHEval80<-c(2.57,2.600,2.615,2.613)
CHEval128<-c(2.57,2.605,2.620,2.62)
SGBMval<-c(2.58,2.60,2.620,2.623)
xdt<-c(0.1,0.05,0.02,0.01)
plot(xdt, CHEval40, type="b", ylim=c(2.5,2.65), xlab="dt", ylab="EPE")
lines(xdt, CHEval80, type="b", col="red")
lines(xdt, CHEval128, type="b", col="blue")
lines(xdt, SGBMval, type="b", col="orange")
legend(x="bottomleft", legend=c("CHE40","CHE80","CHE110","SGBM"), col=c("black", "red", "blue", "orange"), lwd=1,bg="white")

#LSM V0
LSM3<-c(6.215,6.217, 6.232,6.221)
#LSM5<-c(2.679,2.680,2.676)
LSM8<-c(6.215,6.232,6.232,6.220)
SGBMval<-c(6.199,6.200,6.199,6.199)
xdt<-c(0.1,0.05,0.02,0.01)
plot(xdt, LSM3, type="b", ylim=c(6.15,6.25), xlab="dt", ylab="V0")
#lines(xdt, LSM5, type="b", col="red")
lines(xdt, LSM8, type="b", col="blue")
lines(xdt, SGBMval, type="b", col="orange")
legend(x="bottomleft", legend=c("LSM3","LSM8","SGBM"), col=c("black", "blue", "orange"), lwd=1,bg="white")


#LSM EPE
LSM3<-c(2.707,2.741,2.742,2.74)
#LSM5<-c(2.679,2.680,2.676)
LSM8<-c(2.58,2.633,2.647,2.638)
SGBMval<-c(2.58,2.60,2.620,2.623)
xdt<-c(0.1,0.05,0.02,0.01)
plot(xdt, LSM3, type="b", ylim=c(2.55,2.8), xlab="dt", ylab="EPE")
#lines(xdt, LSM5, type="b", col="red")
lines(xdt, LSM8, type="b", col="blue")
lines(xdt, SGBMval, type="b", col="orange")
legend(x="bottomleft", legend=c("LSM3","LSM8","SGBM"), col=c("black", "blue", "orange"), lwd=1,bg="white")






############################################## High bias #######################
TT<-200
dt<-0.05
KQ<-10^3*100
KA<-10^3*100

X101<- matrix(data=NA, nrow=KQ, ncol=(TT+1))
Xa101<- matrix(data=NA, nrow=KA, ncol=(TT+1))

set.seed(2023)

#Under Q
sigma<-0.01
lambda<-0.012

#Generating X under Q 
X101[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  X101[,i]<- r((i-1)*dt,i*dt,X101[,i-1],KQ)
}

#Under P
sigma<-0.006 
lambda<-0.008

#Generating X under P
Xa101[,1]<-fm #Time 0
for (i in 2:(TT+1)){
  Xa101[,i]<- r((i-1)*dt,i*dt,Xa101[,i-1],KA)
}


CHE_bias<-matrix(data=NA, nrow=10, ncol=4)
SGBM_bias<-matrix(data=NA, nrow=10, ncol=4)
LSM_bias<-matrix(data=NA, nrow=10, ncol=4)
set.seed(2023)
  
for (j in 1:10){
  sigma<-0.01
  lambda<-0.012
  
  #Generating X under Q 
  X101[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    X101[,i]<- r((i-1)*dt,i*dt,X101[,i-1],KQ)
  }
  
  #Under P
  sigma<-0.006 
  lambda<-0.008
  
  #Generating X under P
  Xa101[,1]<-fm #Time 0
  for (i in 2:(TT+1)){
    Xa101[,i]<- r((i-1)*dt,i*dt,Xa101[,i-1],KA)
  }
  
  #CHE_bias[j,]<-CheVGLAU(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=20, exercisestep=20)[[1]]
  SGBM_bias[j,]<-SGBM_5(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]]
  #LSM_bias[j,]<-LSM_all(X=X101, XaL=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]]
  print(j)
}


mean(CHE_bias[,2])
mean(CHE_bias[,3])
mean(CHE_bias[,4])

mean(SGBM_bias[,1])+c(-1,1)*1.96/sqrt(20)*sd(SGBM_bias[,1])
mean(SGBM_bias[,2])+c(-1,1)*1.96/sqrt(20)*sd(SGBM_bias[,2])
mean(SGBM_bias[,3])+c(-1,1)*1.96/sqrt(20)*sd(SGBM_bias[,3])
mean(SGBM_bias[,4])+c(-1,1)*1.96/sqrt(20)*sd(SGBM_bias[,4])
sd(SGBMB22[,1])

SGBMB55<-SGBM_bias
#SGBMB22<-SGBM_bias

SGBMB55<-SGBMB55[c(1:5),]
SGBMB22<-SGBMB22[c(1:5),]

SGBMB55[1,]
#plot
windows(10,8)
Xlist<-list("B=2"=SGBMB22[,1], "B=5"=SGBMB55[,1])
stripchart(Xlist, vertical=T, ylim=c(6.198,6.29), pch=1, col="red")
abline(6.200,0)

mean(LSM_bias[,2])
mean(LSM_bias[,3])
mean(LSM_bias[,4])


############### Different bundles of SGBM
SGBM_B5<-SGBM_5(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
SGBM_B5[9]
SGBM_B2<-SGBM_2(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
SGBM_B2[9]
SGBM_B2[1]
SGBM_B1<-SGBM_1(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)[[1]] #MANGLER
SGBM_B1

ASD<-SGBM(X=X101, Xa=Xa101, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)

OPPO<-(Xord1>0.01)
mean(OPPO)
OPP1<- (Xord1<(-0.03))
mean(OPP1)

##################### Check if ass holds
step<-125
XordG<-X101[order(X101[,(step-1)]),step] #Make them into X[,step] into order according to step-1
dsa<-SGBM_B2[[2]][order(X101[,(step-1)]),(step)]
dsa10<-ASD[[2]][,(step)]
#Make X into the bundles
Xord1<- XordG[1:(KQ/2)]
Xord2<- XordG[(KQ/2+1):KQ]
dsa1<- dsa[1:(KQ/2)]
dsa2<- dsa[(KQ/2+1):KQ]
hist(Xord1)
hist(Xord2)

#Order them 
dsa1<- dsa1[order(Xord1)]
dsa2<- dsa2[order(Xord2)]
Xord1<- sort(Xord1)
Xord2<- sort(Xord2)


#Data
DATAX1<-matrix(data=c(rep(1,KQ/2),Xord1,Xord1^2), nrow=3, ncol=KQ/2, byrow = T)
DATAX2<-matrix(data=c(rep(1,KQ/2),Xord2,Xord2^2), nrow=3, ncol=(KQ/2), byrow = T)

dd1<-t(DATAX1)%*%SGBM_B2[[11]][,(step*2)]

dd2<-t(DATAX2)%*%SGBM_B2[[11]][,(step*2+1)]


dsa10<-dsa10[order(X101[,step])]
X101step<-sort(X101[,step])
dd<-dd[order(X101[,step], decreasing = F)]
dsa<-dsa[order(X101[,(step)], decreasing = F)]

#Hole graph
windows(10,8)
plot(Xord1,dsa1, xlim=c(-0.08,0.15), ylim=c(0,50), pch=1, col="blue", ylab="Value", xlab="Interest rate")
lines(Xord1,dd1,col="red", lwd=2)
lines(X101step,dsa10, lty=1, col="black", lwd=2)
points(Xord2,dsa2, col=rgb(1,0,0,alpha=1))
points(Xord1,dsa1, col=rgb(1,0,0,alpha=1))
lines(Xord2,dd2,col="black", lwd=2)
legend(x="bottomleft", legend=c("B1: W-hat","B2: W-hat","B1: Approx","B2: Approx"), col=c("blue", "green", "red", "black"), lwd=c(NA,NA,2,2),
       bg="white", pch=c(1,1,NA,NA) )
legend(x="bottomleft", legend=c("Bundle:2", "Bundle:10"), col=c("red", "black"), lwd=c(NA,2),
       bg="white", pch=c(1,NA) )


#Investigation graph
plot(Xord1,(dsa1-dd1), ylim=c(-0.75,0.9), ylab="Error", xlab="Interest rate")
abline(0,0, lty=2, col="red")
plot(Xord2, (dsa2-dd2),ylim=c(-0.75,1.2), ylab="Error", xlab="Interest rate")
abline(0,0, lty=2, col="red")
mean((dsa1-dd1)<0)
mean(abs(dsa1-dd1)<0.1)

################### Approximation error #####################################

#SGBM
SGBM_app_max<-rep(NA,200)

for (i in 1:200){
  SGBM_app_max[i]<-max(abs(SGBM_app[[11]][,i]-SGBM_app[[2]][,i+1]))
}

plot(SGBM_app_max, xlab="timestep", ylab = "Max error", main="SGBM: max error")

Xord<-X10[order(X10[,step]),]
windows(10,8)
plot(X10[order(X10[,151]),151],SGBM_app[[11]][order(X10[,151]),150], type='l', lwd=2, xlab="Interest rate", ylab="Value", main="t=7,5")
lines(X10[order(X10[,151]),151],SGBM_app[[2]][order(X10[,151]),151], col="red", lty=3, lwd=2)
legend(x="topright", legend=c("SGBM-approx","SGBM"), col=c("black", "red"), lty=c(1,3), lwd=2, bg="white")


SGBM_app<-SGBM_approx(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)






############################# Error on cont. value
#dt=0.05

CHE20<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=20, exercisestep=20)
CHE40<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=40, exercisestep=20)
CHE110<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=110, exercisestep=20)
CHE140<-CheVGLAU(X=X10, Xa=Xa10, KQ=KQ, KA=KA, dt=0.05, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1, kk=5, N=140, exercisestep=20) #Hvad er det præcist som går galt?
CHE102[1]
CHE20[1]
CHE40[1]
CHE102[9]
CHE110[1]
CHE140[1]

mean(X101[,500]<(- 0.15))
mean(X101[,500])
mean(X10[,200])
min(X101)
CHE102NYT[[2]][, 182]
SGBM102[[2]][match(q,X10[,2]), 2]
CHE102[[10]]




plot(40,CHE40[[2]][qq[1], 182], xlim=c(35,160), ylim=c(9.2,9.7))
points(80,CHE102[[2]][qq[1], 182], col="red")
points(110,CHE110[[2]][qq[1], 182], col="blue")
points(X10[qq, 182],CHE140[[2]][qq, 182], col="green")

CHE110[[2]][qq, 2]










###################################### Analysis on LSM on number of poly ##########################################
windows(10,8)

#Difference at time 5.95
step<-119
Xord<-X10[order(X10[,step]),]
Cmord<-SGBM102[[2]][order(X10[,step]),]

#New LSM
LSM102ind<-LSM15(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM102ind3<-LSM8(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM_TJEK<-LSM8(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)
LSM_TJEK1[9]
LSM_TJEK1<-LSM15(X=X10, XaL=Xa10, KQ=KQ, KA=KA, dt=dt, Q=4, EX=6, fm=fm, sigma=0.01, lambda=0.012, kp=1)


LSM102ind3[9]
LSM102ind[9]

#poly(3)
CLSMd<-matrix(data=NA, nrow=KQ, ncol=4)
CLSMd[,1]<-1
CLSMd[,2]<-X10[,step]

for (i in 3:4){
  CLSMd[,i]<-CLSMd[,i-1]*CLSMd[,2]
}
CLSM<-CLSMd%*%LSM102[[10]][,step]
CLSMord<- CLSM[order(X10[,step])]

#poly(8)
dim1<- 9
CLSMd1<-matrix(data=NA, nrow=KQ, ncol=dim1)
CLSMd1[,1]<-1
CLSMd1[,2]<-X10[,step]

for (i in 3:dim1){
  CLSMd1[,i]<-CLSMd1[,i-1]*CLSMd1[,2]
}

CLSM1<-CLSMd1%*%LSM102ind3[[10]][,step]
CLSMord1<- CLSM1[order(X10[,step])]


#poly(15)
dim2<- 16
CLSMd2<-matrix(data=NA, nrow=KQ, ncol=dim2)
CLSMd2[,1]<-1
CLSMd2[,2]<-X10[,step]

for (i in 3:dim2){
  CLSMd2[,i]<-CLSMd2[,i-1]*CLSMd2[,2]
}

CLSM2<-CLSMd2%*%LSM102ind[[10]][,step]
CLSMord2<- CLSM2[order(X10[,step])]

#The plot
windows(10,8)
plot(Xord[,step],Cmord[,step], type="l",lwd=2, main="", ylab="Value", xlab = "Interest rate")
grid(nx=NULL)
lines(Xord[,step],CLSMord1, col="red", lwd=2, lty=4)
lines(Xord[,step],CLSMord, col="green", lwd=2, lty=5)
lines(Xord[,step],CLSMord2, col="orange", lwd=2, lty=3)
legend(x="bottomleft", legend=c("SGBM", "LSM(3)", "LSM(8)", "LSM(15)"), col=c("black", "green", "red",  "orange"), lty=c(1,5,4,3), lwd=2, bg="white")

#Table
LSM102ind[1]
LSM102ind3[1]

LSM102ind[9]
LSM102ind3[9]

#plot to show difference in selected path, which are exercised t=6 (t=4) The third exercise period (The first exercise period)
windows(10,8)

extime<-6
experiod<-3
step<- 121
K1<- S(fm,4,10, sigma=0.01, lambda=0.012) #Always this fixed rate
Xord<-X10[order(X10[,step]),]
IN1<- (U(Xord[,step],extime,10,K1,sigma=0.01, lambda=0.012)>0)
mean(IN1) #IN 46.98%



#poly(3)
CLSMd<-matrix(data=NA, nrow=KQ, ncol=4)
CLSMd[,1]<-1
CLSMd[,2]<-X10[,step]

for (i in 3:4){
  CLSMd[,i]<-CLSMd[,i-1]*CLSMd[,2]
}
CLSM<-CLSMd%*%LSM102[[14]][,experiod]
CLSMord<- CLSM[order(X10[,step])]

#poly(8)
dim1<- 9
CLSMd1<-matrix(data=NA, nrow=KQ, ncol=dim1)
CLSMd1[,1]<-1
CLSMd1[,2]<-X10[,step]

for (i in 3:dim1){
  CLSMd1[,i]<-CLSMd1[,i-1]*CLSMd1[,2]
}

CLSM1<-CLSMd1%*%LSM102ind3[[14]][,experiod]
CLSMord1<- CLSM1[order(X10[,step])]


#poly(15)
dim2<- 16
CLSMd2<-matrix(data=NA, nrow=KQ, ncol=dim2)
CLSMd2[,1]<-1
CLSMd2[,2]<-X10[,step]

for (i in 3:dim2){
  CLSMd2[,i]<-CLSMd2[,i-1]*CLSMd2[,2]
}

CLSM2<-CLSMd2%*%LSM102ind[[14]][,experiod]
CLSMord2<- CLSM2[order(X10[,step])]

plot(Xord[IN1,step],CLSMord2[IN1], ylim=c(0,80), col="orange", ylab="Value", xlab="Interest rate", main="")
points(Xord[IN1,step],CLSMord[IN1], col="green")
points(Xord[IN1,step],CLSMord1[IN1], col="red")
lines(Xord[IN1,step],U(Xord[IN1,step],extime,10,K1,sigma=0.01, lambda=0.012), lwd=2)
legend(x="bottomleft", legend=c("Exercise value", "LSM(3)", "LSM(8)", "LSM(15)"), col=c("black", "green", "red",  "orange"), lty=c(1,1,1,1), lwd=2, bg="white")



#Test
Xordin<-Xord[IN1,step]
mean(CLSMord2[IN1]>U(Xord[IN1,step],6,10,K1,sigma=0.01, lambda=0.012))
max(Xordin[CLSMord2[IN1]<U(Xord[IN1,step],6,10,K1,sigma=0.01, lambda=0.012)])
max(Xordin)
mean(Xordin<0.00206) #ITM-Paths are exercise Exercise 72%








#plot to show difference in selected path, which are exercised (t=4) (The first exercise period)
windows(10,8)

extime<-4
experiod<-1
step<- 81
K1<- S(fm,4,10, sigma=0.01, lambda=0.012) #Always this fixed rate
Xord<-X10[order(X10[,step]),]
IN1<- (U(Xord[,step],extime,10,K1,sigma=0.01, lambda=0.012)>0)
mean(IN1)



#poly(3)
CLSMd<-matrix(data=NA, nrow=KQ, ncol=4)
CLSMd[,1]<-1
CLSMd[,2]<-X10[,step]

for (i in 3:4){
  CLSMd[,i]<-CLSMd[,i-1]*CLSMd[,2]
}
CLSM<-CLSMd%*%LSM102[[14]][,experiod]
CLSMord<- CLSM[order(X10[,step])]

#poly(8)
dim1<- 9
CLSMd1<-matrix(data=NA, nrow=KQ, ncol=dim1)
CLSMd1[,1]<-1
CLSMd1[,2]<-X10[,step]

for (i in 3:dim1){
  CLSMd1[,i]<-CLSMd1[,i-1]*CLSMd1[,2]
}

CLSM1<-CLSMd1%*%LSM102ind3[[14]][,experiod]
CLSMord1<- CLSM1[order(X10[,step])]


#poly(15)
dim2<- 16
CLSMd2<-matrix(data=NA, nrow=KQ, ncol=dim2)
CLSMd2[,1]<-1
CLSMd2[,2]<-X10[,step]

for (i in 3:dim2){
  CLSMd2[,i]<-CLSMd2[,i-1]*CLSMd2[,2]
}

CLSM2<-CLSMd2%*%LSM102ind[[14]][,experiod]
CLSMord2<- CLSM2[order(X10[,step])]

plot(Xord[IN1,step],CLSMord2[IN1], ylim=c(0,80), col="orange", ylab="Value", xlab="Interest rate", main="")
points(Xord[IN1,step],CLSMord[IN1], col="green")
points(Xord[IN1,step],CLSMord1[IN1], col="red")
lines(Xord[IN1,step],U(Xord[IN1,step],extime,10,K1,sigma=0.01, lambda=0.012), lwd=2)
legend(x="bottomleft", legend=c("Exercise value", "LSM(3)", "LSM(8)", "LSM(15)"), col=c("black", "green", "red",  "orange"), lty=c(1,1,1,1), lwd=2, bg="white")



#Test
Xordin<-Xord[IN1,step]
mean(CLSMord2[IN1]>U(Xord[IN1,step],4,10,K1,sigma=0.01, lambda=0.012))
max(Xordin[CLSMord2[IN1]<U(Xord[IN1,step],4,10,K1,sigma=0.01, lambda=0.012)])
max(Xordin)
mean(Xordin<0.00206) #ITM-Paths are exercise Exercise 72%






#plot to show difference in continuation value  (t=3.95) (Where PFE is highest)
windows(10,8)

step<- 80
K1<- S(fm,4,10, sigma=0.01, lambda=0.012)
Xord<-X10[order(X10[,step]),]
IN1<- (U(Xord[,step],6,10,K1,sigma=0.01, lambda=0.012)>0)
mean(IN1) #0.485

#poly(3)
CLSMd<-matrix(data=NA, nrow=KQ, ncol=4)
CLSMd[,1]<-1
CLSMd[,2]<-X10[,step]

for (i in 3:4){
  CLSMd[,i]<-CLSMd[,i-1]*CLSMd[,2]
}
CLSM<-CLSMd%*%LSM102[[10]][,step]
CLSMord<- CLSM[order(X10[,step])]

#poly(8)
dim1<- 9
CLSMd1<-matrix(data=NA, nrow=KQ, ncol=dim1)
CLSMd1[,1]<-1
CLSMd1[,2]<-X10[,step]

for (i in 3:dim1){
  CLSMd1[,i]<-CLSMd1[,i-1]*CLSMd1[,2]
}

CLSM1<-CLSMd1%*%LSM102ind3[[10]][,step]
CLSMord1<- CLSM1[order(X10[,step])]


#poly(15)
dim2<- 16
CLSMd2<-matrix(data=NA, nrow=KQ, ncol=dim2)
CLSMd2[,1]<-1
CLSMd2[,2]<-X10[,step]

for (i in 3:dim2){
  CLSMd2[,i]<-CLSMd2[,i-1]*CLSMd2[,2]
}

CLSM2<-CLSMd2%*%LSM102ind[[10]][,step]
CLSMord2<- CLSM2[order(X10[,step])]

Cmord<-SGBM102[[2]][order(X10[,step]),]


plot(Xord[,step],CLSMord2[], ylim=c(0,80), col="orange", ylab="Value", xlab="Interest rate", main="")
points(Xord[,step],CLSMord[], col="green")
points(Xord[,step],CLSMord1[], col="red")
#lines(Xord[,step],U(Xord[,step],4,10,K1,sigma=0.01, lambda=0.012), lwd=2)
lines(Xord[,step],Cmord[,step])
legend(x="bottomleft", legend=c("LSM(3)", "LSM(8)", "LSM(15)", "SGBM"), col=c("green", "red",  "orange", "black"), lty=c(1,1,1), lwd=2, bg="white")


#ZOOM Ej brug MIGHT BE OF USE
windows(10,8)

plot(Xord[IN1,step],CLSMord2[IN1], ylim=c(3,6), xlim=c(0.002,0.0024))
points(Xord[IN1,step],CLSMord[IN1], col="red")
points(Xord[IN1,step],CLSMord1[IN1], col="green")
lines(Xord[IN1,step],U(Xord[IN1,step],6,10,K1,sigma=0.01, lambda=0.012), lwd=2)



### Table god :=)!!
extime<-6
step<-121
EXX1<-(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)==LSM102ind[[2]][,step])
EXX1<- as.logical(EXX1*(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)!=0))
EXX2<-(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)==LSM102ind3[[2]][,step])
EXX2<- as.logical(EXX2*(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)!=0))
EXX3<-(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)==LSM102[[2]][,step])
EXX3<- as.logical(EXX3*(U(X10[,step],extime,10,K1,sigma=0.01, lambda=0.012)!=0))
mean(EXX3)
mean(EXX2)
mean(EXX1)
mean(EXX1!=EXX2)
mean(EXX1!=EXX3)
mean(EXX3!=EXX2)
EXX3<-(U(X10[,81],40,10,K1,sigma=0.01, lambda=0.012)==LSM102[[2]][,81])


plot(X10[,121],LSM102[[2]][,121])
points(X10[EXX1,121], LSM102ind[[2]][EXX1,121], col="red")

#ZOOM
plot(X10[,121],LSM102[[2]][,121], xlim=c(0.00,0.0025), ylim=c(0,25))
points(X10[,121], LSM102ind[[2]][,121], col="red")
points(X10[,121], LSM102ind3[[2]][,121], col="green")

### End of segment.






