library(collapse)
#source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Swaption/Hull-White Swaption.R")
#Needs Hull-White Swaption.R


#LSM
#Functions

LSM<-function(X,XaL,KQ,KA,dt,Q,EX,fm,sigma,lambda,kp){

  
#Parameters
TT<-(EX+Q)/dt
exercisestep<-1/dt
end<-dt*TT
K<- S(fm,Q,(EX+Q), sigma, lambda)*kp

#Under Q
BL<- matrix(data=NA, nrow=4, ncol=TT)
EQL<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
V<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
bex<- matrix(data=NA, nrow=4, ncol=EX)

#Under P
EAL<-matrix(data=NA, nrow=KA, ncol=(TT+1))



########################### SIMULATION BEGINS #############################################################

startLSM<-Sys.time()
##### LAST PERIOD
time<-TT+1
V[,time]<-U(X[,time],end,end, K,sigma, lambda)
EQL[,time]<-0 #No exposure at end
Y<- X^2
Z<- X*Y

for (j in 1:EX){
  #Intermediate dates
  q<-matrix(data=NA, nrow=exercisestep, ncol=KQ)
  q[exercisestep,]<-V[,time]
  

  for (i in 1:(exercisestep-1)){
    q[(exercisestep-i),]<- q[(exercisestep+1-i),]*exp(-0.5*(X[,(time-1)]+X[,(time)])*dt)
    
    BL[,time-1]<- flm(X=cbind(1,X[,time-1],Y[,time-1],Z[,time-1]),y=q[exercisestep-i,])
    EQL[,time-1]<- BL[1,time-1]+BL[2,time-1]*X[,time-1]+BL[3,time-1]*Y[,time-1]+BL[4,time-1]*Z[,time-1]
    time<-time-1}
  
  #Exercising
  time<- time -1 #Since notation is plus 1
  rtime<- (time-1)*dt
  IN<- (U(X[,time],rtime,end, K,sigma, lambda)>0)
  
  XINN<-X[IN,time]
  XINN2<-X[IN,time+1]
  VINN<- q[1,IN]*exp(-0.5*(XINN+XINN2)*dt)
  VEX<-q[1,]*exp(-0.5*(X[,time]+X[,time+1])*dt)
  
  
  YINN<- XINN^2
  ZINN <- YINN*XINN
  

  bex[,EX-j+1]<-flm(X=cbind(1,XINN,YINN,ZINN),y=VINN)
  
  ###
  
  CASHH<- bex[1,EX-j+1]+bex[2,EX-j+1]*XINN+bex[3,EX-j+1]*YINN+bex[4,EX-j+1]*ZINN #Cashflow for ITM options
  NEXX <- CASHH>U(XINN,rtime,end, K,sigma, lambda) #Not exercising
  EXX <- CASHH<U(XINN,rtime,end, K,sigma, lambda)
  qIN<-q[1,IN]
  VINN1<- U(XINN,rtime,end, K,sigma, lambda) #"ALL EXERCISING"
  VINN1[NEXX]<- exp(-0.5*(XINN[NEXX]+XINN2[NEXX])*dt)*qIN[NEXX] #Those who are not exercised
  
  # 
  V[,time]<- exp(-0.5*(X[,time]+X[,time+1])*dt)*q[1,] #All
  V[IN,time]<- VINN1
  
  #EQL
  BL[,time]<-flm(X=cbind(1,X[,time],Y[,time],Z[,time]),y=VEX)
  EQL[,time]<-BL[1,time]+BL[2,time]*X[,time]+BL[3,time]*Y[,time]+BL[4,time]*Z[,time] #"NON EXERCISING"
  EQLIN <- EQL[IN,time]
  EQLIN[EXX]<- 0
  EQL[IN,time]<-EQLIN
  EQL[EQL[,time]==0,]<-0
}



#Without exercising
for (p in 1:Q){
  q<-matrix(data=NA, nrow=exercisestep, ncol=KQ)
  q[exercisestep,]<-V[,time]
  
  for (i in 1:(exercisestep-1)){
    q[(exercisestep-i),]<- q[(exercisestep+1-i),]*exp(-0.5*(X[,(time-1)]+X[,(time)])*dt)
    
    BL[,time-1]<- flm(X=cbind(1,X[,time-1],Y[,time-1],Z[,time-1]),y=q[exercisestep-i,])
    EQL[,time-1]<- BL[1,time-1]+BL[2,time-1]*X[,time-1]+BL[3,time-1]*Y[,time-1]+BL[4,time-1]*Z[,time-1]
    time<-time-1
  }
  
  time<- time-1
  
  V[,time]<-exp(-0.5*(X[,time]+X[,time+1])*dt)*q[1,]
  BL[,time]<- flm(X=cbind(1,X[,time],Y[,time],Z[,time]),y=V[,time])
  EQL[,time]<- BL[1,time]+BL[2,time]*X[,time]+BL[3,time]*Y[,time]+BL[4,time]*Z[,time]
}

EQL[,time]<- BL[1,time]






##### REAL WORLD

startal<-Sys.time()
YaL<-XaL^2
ZaL<-YaL*XaL


##### LAST PERIOD
time<-TT+1
EAL[,time]<-0

for (j in 1:EX){
  #Intermediate dates
  
  for (i in 1:(exercisestep-1)){
    EAL[,time-1]<- BL[1,time-1]+BL[2,time-1]*XaL[,time-1]+BL[3,time-1]*YaL[,time-1]+BL[4,time-1]*ZaL[,time-1]
    time<-time-1
  }
  
  #Exercising
  time<- time -1 
  rtime<- (time-1)*dt
  
  
  IN<- (U(XaL[,time],rtime,end, K, sigma, lambda)>0)
  XINN<-XaL[IN,time]
  YINN<- XINN^2
  ZINN <- YINN*XINN
  ###
  
  Co<- BL[1,time]+BL[2,time]*XaL[,time]+BL[3,time]*YaL[,time]+BL[4,time]*ZaL[,time]
  EAL[,time]<-Co
  EALIN <- EAL[IN,time]
  CASHH<-bex[1,EX-j+1]+bex[2,EX-j+1]*XINN+bex[3,EX-j+1]*YINN+bex[4,EX-j+1]*ZINN
  EXX <- CASHH<U(XINN,rtime,end, K,sigma, lambda) #Exercising
  EALIN[EXX]<-0
  EAL[IN,time]<-EALIN
  EAL[EAL[,time]==0,]<-0
  
}

#Without exercising
for (p in 1:(Q/dt-1)){
  EAL[,time-1]<- BL[1,time-1]+BL[2,time-1]*XaL[,time-1]+BL[3,time-1]*YaL[,time-1]+BL[4,time-1]*ZaL[,time-1]
  time<-time-1
}

EAL[,time-1]<-mean(V[,time-1])

slutLSM<-Sys.time()



########################## Estimation ########################

#Under P
PFEL<-rep(NA,TT+1)

for (i in 1:(TT+1)){
  PFEL[i]<-quantile(EAL[,i], p=0.99)
}

mEAL<-rep(NA,TT+1)
for (i in 1:(TT+1)){
  mEAL[i]<-mean(EAL[,i])
}


#Under Q

#Discount EQ
DEQL<-matrix(data=NA, nrow=KQ, ncol=TT+1)
mDEQL<-rep(NA,TT+1)
mEQL<-rep(NA,TT+1)
PFEql<-rep(NA,TT+1)

AQL<- 0
i<-1
DEQL[,i]<- exp(AQL)*EQL[,i]
mDEQL[i]<-mean(DEQL[,i])
PFEql[i]<- quantile(EQL[,i],p=0.99)
mEQL[i]<-mean(EQL[,i])


for (i in 2:(TT+1)){
  AQL<- AQL- 0.5*(X[,i-1]+X[,i])*dt
  DEQL[,i]<- exp(AQL)*EQL[,i]
  mDEQL[i]<-mean(DEQL[,i])
  mEQL[i]<-mean(EQL[,i])
  PFEql[i]<- quantile(EQL[,i],p=0.99)
}

#CVA
dPS<-function(t){
  1-exp(-0.02*t)
}
point<-seq(0,end,dt)


CVAL<-rep(NA,TT)


for (i in 1:(TT)){
  CVAL[i]<- mDEQL[i+1]*(dPS(point[i+1])-dPS(point[i]))
}








### Creating output

output<-list()

output[[1]]<-c(mean(V[,1]), max(PFEL), mean(Re(mEAL[-1])),100*sum(CVAL))
output[[2]]<-V
output[[3]]<-EQL
output[[4]]<-EAL
output[[5]]<-PFEL
output[[6]]<-PFEql
output[[7]]<-mEAL
output[[8]]<- mEQL
output[[9]]<- slutLSM-startLSM
output[[10]]<-BL
output[[11]]<- X
output[[12]]<- Y
output[[13]]<- Z
output[[14]]<-bex

names(output)<-c("V0, MPFE, EPE, 100CVA", "V", "EQ", "EA", "PFE", "PFEql",  "mEA", "mEQ", "Tid", "BL", "X", "Y", "Z", "bex")

return(output)
}
