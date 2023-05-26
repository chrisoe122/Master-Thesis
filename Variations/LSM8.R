library(dplyr)
library(forestmangr)
library(data.table)
library(geometry)
library(plyr)
library(RcppEigen)
library(collapse)



#LSM
#Functions


source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Swaption/Hull-White Swaption.R")


LSM8<-function(X,XaL,KQ,KA,dt,Q,EX,fm,sigma,lambda,kp){
  
  
  #Parameters
  TT<-(EX+Q)/dt
  exercisestep<-1/dt ############### HARDCODED
  end<-dt*TT
  K<- S(fm,Q,(EX+Q), sigma, lambda)*kp
  
  #Under Q
  BL<- matrix(data=NA, nrow=9, ncol=TT)
  EQL<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
  V<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
  bex<- matrix(data=NA, nrow=9, ncol=EX)
  
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
  Z1<- X*Z
  Z2<- X*Z1
  Z3<- X*Z2
  Z4<- X*Z3
  Z5<- X*Z4
  
  for (j in 1:EX){
    #Intermediate dates
    q<-matrix(data=NA, nrow=exercisestep, ncol=KQ)
    q[exercisestep,]<-V[,time]
    
    
    for (i in 1:(exercisestep-1)){
      q[(exercisestep-i),]<- q[(exercisestep+1-i),]*exp(-0.5*(X[,(time-1)]+X[,(time)])*dt)
      
      BL[,time-1]<- flm(X=cbind(1,X[,time-1],Y[,time-1],Z[,time-1],Z1[,time-1],Z2[,time-1],Z3[,time-1],Z4[,time-1],Z5[,time-1]),y=q[exercisestep-i,])
      EQL[,time-1]<- BL[1,time-1]+BL[2,time-1]*X[,time-1]+BL[3,time-1]*Y[,time-1]+BL[4,time-1]*Z[,time-1]+
        BL[5,time-1]*Z1[,time-1]+BL[6,time-1]*Z2[,time-1]+BL[7,time-1]*Z3[,time-1]+BL[8,time-1]*Z4[,time-1]+BL[9,time-1]*Z5[,time-1]
        
      time<-time-1}
    
    #Exercising
    time<- time -1 #Da notation her er plus 1
    rtime<- (time-1)*dt
    IN<- (U(X[,time],rtime,end, K,sigma, lambda)>0)
    
    XINN<-X[IN,time]
    XINN2<-X[IN,time+1]
    VINN<- q[1,IN]*exp(-0.5*(XINN+XINN2)*dt)
    VEX<-q[1,]*exp(-0.5*(X[,time]+X[,time+1])*dt)
    
    
    YINN<- XINN^2
    ZINN <- YINN*XINN
    ZINN1 <- ZINN*XINN
    ZINN2 <- ZINN1*XINN
    ZINN3 <- ZINN2*XINN
    ZINN4 <- ZINN3*XINN
    ZINN5 <- ZINN4*XINN
    
    
    bex[,EX-j+1]<-flm(X=cbind(1,XINN,YINN,ZINN,ZINN1,ZINN2,ZINN3,ZINN4,ZINN5),y=VINN)
    
    ###
    
    CASHH<- bex[1,EX-j+1]+bex[2,EX-j+1]*XINN+bex[3,EX-j+1]*YINN+bex[4,EX-j+1]*ZINN+
      bex[5,EX-j+1]*ZINN1+bex[6,EX-j+1]*ZINN2+bex[7,EX-j+1]*ZINN3+bex[8,EX-j+1]*ZINN4+bex[9,EX-j+1]*ZINN5#Cashflow for ITM options
    NEXX <- CASHH>U(XINN,rtime,end, K,sigma, lambda) #Not exercising
    EXX <- CASHH<U(XINN,rtime,end, K,sigma, lambda)
    qIN<-q[1,IN]
    VINN1<- U(XINN,rtime,end, K,sigma, lambda) #"ALL EXERCISING"
    VINN1[NEXX]<- exp(-0.5*(XINN[NEXX]+XINN2[NEXX])*dt)*qIN[NEXX] #Those who are not exercised
    #VINN1[NEXX]<- CASHH[NEXX] #Those who are not exercised
    
    # 
    V[,time]<- exp(-0.5*(X[,time]+X[,time+1])*dt)*q[1,] #All
    V[IN,time]<- VINN1
    
    #EQL
    BL[,time]<-flm(X=cbind(1,X[,time],Y[,time],Z[,time],Z1[,time],Z2[,time],Z3[,time],Z4[,time],Z5[,time]),y=VEX)
    EQL[,time]<-BL[1,time]+BL[2,time]*X[,time]+BL[3,time]*Y[,time]+BL[4,time]*Z[,time]+
      BL[5,time]*Z1[,time]+BL[6,time]*Z2[,time]+BL[7,time]*Z3[,time]+BL[8,time]*Z4[,time]+BL[9,time]*Z5[,time] #"NON EXERCISING"
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
      
      BL[,time-1]<- flm(X=cbind(1,X[,time-1],Y[,time-1],Z[,time-1],Z1[,time-1],Z2[,time-1],Z3[,time-1],Z4[,time-1],Z5[,time-1]),y=q[exercisestep-i,])
      EQL[,time-1]<- BL[1,time-1]+BL[2,time-1]*X[,time-1]+BL[3,time-1]*Y[,time-1]+BL[4,time-1]*Z[,time-1]+
        BL[5,time-1]*Z1[,time-1]+BL[6,time-1]*Z2[,time-1]+BL[7,time-1]*Z3[,time-1]+BL[8,time-1]*Z4[,time-1]+BL[9,time-1]*Z5[,time-1]
      time<-time-1
    }
    
    time<- time-1
    
    V[,time]<-exp(-0.5*(X[,time]+X[,time+1])*dt)*q[1,]
    BL[,time]<-flm(X=cbind(1,X[,time],Y[,time],Z[,time],Z1[,time],Z2[,time],Z3[,time],Z4[,time],Z5[,time]),y=VEX)
    EQL[,time]<-BL[1,time]+BL[2,time]*X[,time]+BL[3,time]*Y[,time]+BL[4,time]*Z[,time]+
      BL[5,time]*Z1[,time]+BL[6,time]*Z2[,time]+BL[7,time]*Z3[,time]+BL[8,time]*Z4[,time]+BL[9,time]*Z5[,time] #"NON EXERCISING"
  }
  
  EQL[,time]<- BL[1,time]
  
  
  
  
  
  
  ##### REAL WORLD
  
  startal<-Sys.time()
  YaL<-XaL^2
  ZaL<-YaL*XaL
  ZaL1<- ZaL*XaL
  ZaL2<- ZaL1*XaL
  ZaL3<- ZaL2*XaL
  ZaL4<- ZaL3*XaL
  ZaL5<- ZaL4*XaL
  
  
  ##### LAST PERIOD
  time<-TT+1
  EAL[,time]<-0
  
  for (j in 1:EX){
    #Intermediate dates
    
    for (i in 1:(exercisestep-1)){
      EAL[,time-1]<- BL[1,time-1]+BL[2,time-1]*XaL[,time-1]+BL[3,time-1]*YaL[,time-1]+BL[4,time-1]*ZaL[,time-1]+
        BL[5,time-1]*ZaL1[,time-1]+BL[6,time-1]*ZaL2[,time-1]+BL[7,time-1]*ZaL3[,time-1]+BL[8,time-1]*ZaL4[,time-1]+BL[9,time-1]*ZaL5[,time-1]
      time<-time-1
    }
    
    #Exercising
    time<- time -1 #Da notation her er plus 1
    rtime<- (time-1)*dt
    
    
    IN<- (U(XaL[,time],rtime,end, K, sigma, lambda)>0)
    XINN<-XaL[IN,time]
    YINN<- XINN^2
    ZINN <- YINN*XINN
    ZINN1 <- ZINN*XINN
    ZINN2 <- ZINN1*XINN
    ZINN3 <- ZINN2*XINN
    ZINN4 <- ZINN3*XINN
    ZINN5 <- ZINN4*XINN
    ###
    
    Co<- BL[1,time]+BL[2,time]*XaL[,time]+BL[3,time]*YaL[,time]+BL[4,time]*ZaL[,time]+
      BL[5,time]*ZaL1[,time]+BL[6,time]*ZaL2[,time]+BL[7,time]*ZaL3[,time]+BL[8,time]*ZaL4[,time]+BL[9,time]*ZaL5[,time]
    EAL[,time]<-Co
    EALIN <- EAL[IN,time]
    CASHH<-bex[1,EX-j+1]+bex[2,EX-j+1]*XINN+bex[3,EX-j+1]*YINN+bex[4,EX-j+1]*ZINN+
      bex[5,EX-j+1]*ZINN1+bex[6,EX-j+1]*ZINN2+bex[7,EX-j+1]*ZINN3+bex[8,EX-j+1]*ZINN4+bex[9,EX-j+1]*ZINN5
    EXX <- CASHH<U(XINN,rtime,end, K,sigma, lambda) #Exercising
    EALIN[EXX]<-0
    EAL[IN,time]<-EALIN
    EAL[EAL[,time]==0,]<-0
    
  }
  
  #Without exercising
  for (p in 1:(Q/dt-1)){
    EAL[,time-1]<- BL[1,time-1]+BL[2,time-1]*XaL[,time-1]+BL[3,time-1]*YaL[,time-1]+BL[4,time-1]*ZaL[,time-1]+
      BL[5,time-1]*ZaL1[,time-1]+BL[6,time-1]*ZaL2[,time-1]+BL[7,time-1]*ZaL3[,time-1]+BL[8,time-1]*ZaL4[,time-1]+BL[9,time-1]*ZaL5[,time-1]
    time<-time-1
  }
  
  EAL[,time-1]<-mean(V[,time-1]) #Eksponering i første periode, må være værdien
  
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
  
  
  for (i in 1:(TT)){ #Skal 0 være med?
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
