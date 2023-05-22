library(dplyr)
library(forestmangr)
library(data.table)
library(geometry)
library(plyr)
library(RcppEigen)
library(collapse)



#SGBM

source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Swaption/Hull-White Swaption.R")






SGBM<-function(X,Xa,KQ,KA,dt,Q,EX,fm,sigma1,lambda1,kp){

  ### Function to create charistarics function and its derivative
  ###Charastic function
  Atilde1<-function(tau){
    sigma^2/(2*lambda^3)*(lambda*tau-2*(1-exp(-lambda*tau))+1/2*(1-exp(-2*lambda*tau)))
  } #PAPER have changed sigma og lambda
  
  
  Btilde1<-function(tau){
    -1/lambda*(1-exp(-lambda*tau))
  }
  
  thetaint<-function(t,TT){
    fm*(TT-t)+sigma^2/(4*lambda^3)*(2*lambda*(TT-t)+4*(exp(-lambda*TT)-exp(-lambda*t))-(exp(-2*lambda*TT)-exp(-2*lambda*t)))
  } #PAPER (QIAN) have forgotten a  2 in front of lambda.
  
  
  theta<-function(t){
    fm+sigma^2/(2*lambda^2)*(1-exp(-lambda*t))^2
  }
  
  dChf1<-function(t,TT,r){ #Chr. fnc
    exp(-thetaint(t,TT)+Atilde1(TT-t)+Btilde1(TT-t)*(r-theta(t)))
  } #VALIDERET
  
  difChf1<-function(t,TT,r){# Diff chr. fnc.
    dChf1(t,TT,r)*(theta(TT)-sigma^2/(2*lambda^2)*(1-exp(-lambda*(TT-t)))^2+
                     exp(-lambda*(TT-t))*(r-theta(t)))
  } #VALIDERET
  
  
  dif2Chf1<-function(t,TT,r){# Double diff chr. fnc
    dChf1(t,TT,r)*sigma^2/(2*lambda)*(1-exp(-2*lambda*(TT-t)))+
      dChf1(t,TT,r)*(theta(TT)-sigma^2/(2*lambda^2)*(1-exp(-lambda*(TT-t)))^2+exp(-lambda*(TT-t))*(r-theta(t)))^2
  } #VALIDERET
  
  
  
### Parameters
sigma<-sigma1
lambda<-lambda1
K<- S(fm,Q,Q+EX, sigma1, lambda1)*kp #ATM strike price times a procent
TT<-(Q+EX)/dt #Time-periods
end<-TT*dt #End time

#Under Q
bmny<-matrix(data=NA, nrow=3, ncol=((TT-1)*10+1)) #col er tidsgrupperinger med 10 bundles i hver, så TT*10. row er antal parameter
Cm<-matrix(data=NA, nrow=KQ, ncol=TT) 
EQ<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
Dm<-matrix(data=NA, nrow=9, ncol=(TT)) #Upper bound for hver bundle

#Under P
Cma<-matrix(data=NA, nrow=KA, ncol=TT) 
EA<-matrix(data=NA, nrow=KA, ncol=(TT+1))





########################### SIMULATION BEGINS #############################################################

startSG<- Sys.time()
##### SGBM

# LAST PERIOD
time<-TT+1
Vny<-matrix(data=NA, nrow=KQ, ncol=(TT+1))
EQ[,time]<-0

#Generate the last V
Vny[,(TT+1)]<-U(X[,time],end,end, K, sigma1, lambda1)

for( p in 1:EX){
  
  ###Intermediates dates
  #dsad<-Sys.time()
  for (j in 1:(1/dt-1)){ #ÆNDRER
    
    #Calculations
    df1<-data.table(Vny[,time],X[,time], X[,time]^2, X[,time-1])
    df1$quan<-ntile(df1[[4]], 10)#Number of groups
    CoefSG<-df1[,as.list(flm(X=cbind(1,V2,V3),y=V1)),by="quan"]
    bmny[,((time-3)*10+2):((time-2)*10+1)]<- t(CoefSG[order(quan),-1])
    df1<-inner_join(df1,CoefSG, by="quan")
    
    læ<-dChf1((time-2)*dt,(time-1)*dt,df1[[4]])
    læ1<-difChf1((time-2)*dt,(time-1)*dt,df1[[4]])
    læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,df1[[4]])
    dChrdata<-rbind(læ,læ1,læ2)
    kle<- rbind(df1[[6]], df1[[7]], df1[[8]])
    
    #Update values
    Cm[,(time-1)]<- dot(kle,dChrdata)
    Dm[,(time-1)]<-quantile(df1[[4]],probs = seq(0.1,0.9,0.1))
    EQ[,time-1]<- Re(Cm[,time-1])
    Vny[,time-1]<- Re(Cm[,time-1])
    time<-time-1
  }
  ### Period with exercise possiblility
  
  #Calculations
  df1<-data.table(Vny[,time],X[,time], X[,time]^2, X[,time-1])
  df1$quan<-ntile(df1[[4]], 10)#Number of groups
  CoefSG<-df1[,as.list(flm(X=cbind(1,V2,V3),y=V1)),by="quan"]
  bmny[,((time-3)*10+2):((time-2)*10+1)]<- t(CoefSG[order(quan),-1])
  df1<-inner_join(df1,CoefSG, by="quan")
  
  læ<-dChf1((time-2)*dt,(time-1)*dt,df1[[4]])
  læ1<-difChf1((time-2)*dt,(time-1)*dt,df1[[4]])
  læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,df1[[4]])
  dChrdata<-rbind(læ,læ1,læ2)
  kle<- rbind(df1[[6]], df1[[7]], df1[[8]])
  
  #Update values
  Cm[,(time-1)]<- dot(kle,dChrdata)
  Dm[,(time-1)]<-quantile(df1[[4]],probs = seq(0.1,0.9,0.1))
  Vny[,time-1]<- pmax(Re(Cm[,time-1]),U(X[,time-1],(time-2)*dt,end,K,sigma1, lambda1))
  EQ[,time-1]<-Re(Cm[,time-1]) 
  EQ[(U(X[,time-1],(time-2)*dt,end,K,sigma1, lambda1)==Vny[,time-1]),]<-0
  time<- time-1
}

###Intermediates dates
for (j in 1:(Q/dt-1)){
  
  #Calculations
  df1<-data.table(Vny[,time],X[,time], X[,time]^2, X[,time-1])
  df1$quan<-ntile(df1[[4]], 10)#Number of groups
  CoefSG<-df1[,as.list(flm(X=cbind(1,V2,V3),y=V1)),by="quan"]
  bmny[,((time-3)*10+2):((time-2)*10+1)]<- t(CoefSG[order(quan),-1])
  df1<-inner_join(df1,CoefSG, by="quan")
  
  læ<-dChf1((time-2)*dt,(time-1)*dt,df1[[4]])
  læ1<-difChf1((time-2)*dt,(time-1)*dt,df1[[4]])
  læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,df1[[4]])
  dChrdata<-rbind(læ,læ1,læ2)
  kle<- rbind(df1[[6]], df1[[7]], df1[[8]])
  
  #Update values
  Cm[,(time-1)]<- dot(kle,dChrdata)
  Dm[,(time-1)]<-quantile(df1[[4]],probs = seq(0.1,0.9,0.1))
  EQ[,time-1]<- Re(Cm[,time-1])
  Vny[,time-1]<- Re(Cm[,time-1])
  time<-time-1
}

#First period 

df1<-data.table(Vny[,time],X[,time], X[,time]^2, X[,time-1])
CoefSG<-df1[,as.list(flm(X=cbind(1,V2,V3),y=V1))]
bmny[,1]<- t(CoefSG[1,])
læ<-dChf1((time-2)*dt,(time-1)*dt,fm)
læ1<-difChf1((time-2)*dt,(time-1)*dt,fm)
læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,fm)
dChrdata<-rbind(læ,læ1,læ2)

Cm[,(time-1)]<- t(dChrdata)%*%bmny[,1]
EQ[,time-1]<- Re(Cm[,time-1])
Vny[,time-1]<- Re(Cm[,time-1])




############# REAL WORLD  #############



####### LAST PERIOD
time<-TT+1

#Generate the last EA
EA[,TT+1]<-0

qlabels<-c("1","2","3","4","5", "6","7","8","9","10")

for( p in 1:EX){
  
  ###Intermediates dates
  for (j in 1:(1/dt-1)){ #ÆNDRER
    #print((time-1)*dt)
    qbreak<-c(-1,Dm[,(time-1)],1)
    dfa<-data.table(Xa[,time-1]) #Tjek at index passer
    setDT(dfa)[,qlabels:=cut(V1,qbreak, labels=qlabels)]
    dfb<-data.table(t(bmny[,((time-3)*10+2):((time-2)*10+1)]),qlabels)
    dfc<-inner_join(dfa,dfb,by="qlabels")
    
    læ<-dChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
    læ1<-difChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
    læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,dfc[[1]])
    dChrdata<-rbind(læ,læ1,læ2)
    kle<- rbind(dfc[[3]], dfc[[4]], dfc[[5]])
    
    EA[,time-1]<-  dot(kle,dChrdata)
    time<-time-1
  }
  
  ### Period with exercise possiblility
  qbreak<-c(-1,Dm[,(time-1)],1)
  dfa<-data.table(Xa[,time-1])
  setDT(dfa)[,qlabels:=cut(V1,qbreak, labels=qlabels)]
  dfb<-data.table(t(bmny[,((time-3)*10+2):((time-2)*10+1)]),qlabels)
  dfc<-inner_join(dfa,dfb,by="qlabels")
  
  læ<-dChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  læ1<-difChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  dChrdata<-rbind(læ,læ1,læ2)
  kle<- rbind(dfc[[3]], dfc[[4]], dfc[[5]])
  
  Cma[,time-1]<-dot(kle,dChrdata)
  c<-Re(Cma[,time-1])>U(Xa[,time-1],(time-2)*dt,end,K, sigma1, lambda1)
  EA[c,time-1]<-Re(Cma[c,time-1]) 
  d<-Re(Cma[,time-1])<=U(Xa[,time-1],(time-2)*dt,end,K, sigma1, lambda1)
  EA[d,]<-0
  
  
  
  
  time<- time-1
}



for (q in 1:(Q/dt-1)){
  
  ###Intermediates dates
  qbreak<-c(-1,Dm[,(time-1)],1)
  dfa<-data.table(Xa[,time-1])
  setDT(dfa)[,qlabels:=cut(V1,qbreak, labels=qlabels)]
  dfb<-data.table(t(bmny[,((time-3)*10+2):((time-2)*10+1)]),qlabels)
  dfc<-inner_join(dfa,dfb,by="qlabels")
  
  læ<-dChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  læ1<-difChf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,dfc[[1]])
  dChrdata<-rbind(læ,læ1,læ2)
  kle<- rbind(dfc[[3]], dfc[[4]], dfc[[5]])
  
  EA[,time-1]<-  dot(kle,dChrdata)
  time<-time-1
}

#Sidste period

læ<-dChf1((time-2)*dt,(time-1)*dt,fm)
læ1<-difChf1((time-2)*dt,(time-1)*dt,fm)
læ2<-dif2Chf1((time-2)*dt,(time-1)*dt,fm)
dChrdata<-rbind(læ,læ1,læ2)

Cm[,(time-1)]<- t(dChrdata)%*%bmny[,1]  
EA[,time-1]<-t(dChrdata)%*%bmny[,1]
slutSG<-Sys.time()




########################## Estimation ########################

#Under P
PFE<-rep(NA,TT+1)

for (i in 1:(TT+1)){
  PFE[i]<-quantile(Re(EA[,i]), p=0.99)
}

mEA<- rep(NA,TT+1)

for (i in 1:(TT+1)){
  mEA[i]<-mean(EA[,i])
}



#Under Q

#Discount EQ
DEQ<-matrix(data=NA, nrow=KQ, ncol=TT+1)
mEQ<-rep(NA,TT+1)
mDEQ<-rep(NA,TT+1)
PFEq<-rep(NA,TT+1)

AQ<- 0
i<-1
DEQ[,i]<- exp(AQ)*EQ[,i]
mEQ[i]<-mean(EQ[,i])
mDEQ[i]<-mean(DEQ[,i])
PFEq[i]<-quantile(EQ[,i],p=0.99)

for (i in 2:(TT+1)){
  AQ<- AQ- 0.5*(X[,i-1]+X[,i])*dt
  DEQ[,i]<- exp(AQ)*EQ[,i]
  mEQ[i]<-mean(EQ[,i])
  mDEQ[i]<-mean(DEQ[,i])
  PFEq[i]<-quantile(EQ[,i],p=0.99)
}



#CVA
dPS<-function(t){
  1-exp(-0.02*t)
}
point<-seq(0,end,dt)


CVA<-rep(NA,TT)

for (i in 1:(TT)){
  CVA[i]<- mDEQ[i+1]*(dPS(point[i+1])-dPS(point[i]))
}




##### Creating output  ######

output<-list()

output[[1]]<-c(mean(Vny[,1]), max(PFE), mean(Re(mEA[-1])),100*sum(CVA))
output[[2]]<-Vny
output[[3]]<-EQ
output[[4]]<-EA
output[[5]]<-PFE
output[[6]]<- PFEq
output[[7]]<-mEA
output[[8]]<- mEQ
output[[9]]<- slutSG-startSG
output[[10]]<-Cm

names(output)<-c("V0, MPFE, EPE, 100CVA", "V", "EQ", "EA", "PFE", "PFEq", "mEA", "mEQ", "Tid", "C")



return(output)
}
