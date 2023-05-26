#source("C:/Users/Chris/Documents/Mat-øk/5. år/Speciale/R-kode/Swaption/Hull-White Swaption.R")
#Needs Hull-White Swaption.R

CheVGLAU<-function(X,XaChe,KQ,KA,dt,Q,EX,fm,sigma,lambda,kp, kk, N, exercisestep){
  
  sigma<-sigma
  lambda<-lambda
  
  ################ FUNCTIONS
  
  tau<-function(xned,xop,z){
    xop+0.5*(xned-xop)*(1-z)
  }
  
  mun<-function(t,x){
    thetache(t)-lambda*x
  }
  
  tauinv<-function(y,xop1=xop,xned1=xned){
    1-(xop1-y)*2/(xop1-xned1)
  }
  
  thetache<-function(t){
    lambda*fm+sigma^2/(2*lambda)*(1-exp(-2*lambda*t))
  }
  
  
  meanpj<-function(N,mu,sigma){
    #Parameters
    T0<- 1
    T1<- -1
    
    #Y variable
    meany<- tauinv(mu) #Expectation of tauinv 
    vary<- (2/(xop-xned))^2*sigma^2
    
    result<-rep(NA,N+1)
    
    ####mu0
    result[1]<-pnorm(1,mean=meany, sd=sqrt(vary))-pnorm(-1,mean=meany, sd=sqrt(vary))
    
    ###mu1
    result[2]<-meany*result[1]-vary*(dnorm(1,mean=meany, sd=sqrt(vary))-dnorm(-1,mean=meany, sd=sqrt(vary)))
    
    ####mu2
    i<-2
    result[i+1]<- 2*meany*result[i]-2*vary*(dnorm(1,mean=meany, sd=sqrt(vary))-dnorm(-1,mean=meany, sd=sqrt(vary))*T1-
                                              2*(i-1)*((0.5*result[1]*((i-1)%%2==1))))-result[i-1]
    
    #Update T
    T1n<- 2*T1*(-1)-T0
    T0<-T1
    T1<-T1n
    
    
    for (i in 3:N){
      result[i+1]<- 2*meany*result[i]-2*vary*(dnorm(1,mean=meany, sd=sqrt(vary))-dnorm(-1,mean=meany, sd=sqrt(vary))*T1-
                                                2*(i-1)*(0.5*result[1]*((i-1)%%2==1)+sum(result[2:(i-1)]*indi((1:(i-2)),i-1))))-result[i-1]
      T1n<- 2*T1*(-1)-T0
      T0<-T1
      T1<-T1n
      
    }
    return(result)
  }
  
  Gammat<-function(N, mu, sigma){ #Function to create Gamma
    hold<-matrix(data=NA, nrow=length(mu), ncol=(N+1)) #col=p_j, row=x_j
    
    for (i in 1:length(mu)){
      hold[i,]<-meanpj(N,mu[i],sigma)
    }
    return(hold)
  }
  
  indi<-function(j,n){ # Moduls function. Used in the function 'meanpj'
    ((j+n)%%2==1)
  }
  
  alpha<-function(t){ 
    0.01+sigma^2/(2*lambda^2)*(1-exp(-lambda*t))^2 #Comes from Brigo-Interest rate s. 73
  }
  
  
  
  
  
  ############### Chebyshev functions
  #Under Q-measure
  Tj<-function(j,zk){
    cos(j*acos(zk))
  }
  
  CCV<-function(Nk1,time){ ####Depends on z
    NN<-Nk1+1
    cc<-rep(NA,NN)
    
    cc[1]<- 1/(Nk1)*((0.5*(VCh[1,time]*Tj(0,z[1])+VCh[NN,time]*Tj(0,z[NN])))+
                       sum(VCh[-c(1,NN),time]*Tj(0,z[-c(1,NN)])))
    
    for (i in 2:(Nk1)){
      cc[i]<-2*1/(Nk1)*((0.5*(VCh[1,time]*Tj((i-1),z[1])+VCh[NN,time]*Tj((i-1),z[NN])))+
                          sum(VCh[-c(1,NN),time]*Tj((i-1),z[-c(1,NN)])))
    }
    cc[NN]<-1/(Nk1)*((0.5*(VCh[1,time]*Tj(Nk1,z[1])+VCh[NN,time]*Tj(Nk1,z[NN])))+
                       sum(VCh[-c(1,NN),time]*Tj(Nk1,z[-c(1,NN)])))
    return(cc)
  }
  
  IV<-function(x,Nk1,time){
    # Map x into [-1, 1]
    xx1 <- (2*x - (xop+xned))/(xop-xned)
    
    C<-CCV(Nk1,time)
    len<-Nk1+1
    
    #Evaluate the Chebyshev polynomial
    b2<-rep(0,length(x))
    b1<-rep(0,length(x))
    
    for (i in 1:(len-1)){
      b0<- C[len-(i-1)]+2*xx1*b1-b2
      b2<-b1
      b1<-b0
    }
    f<-C[1]+b1*xx1-b2
    
    return(f)
  }
  
  #Under P-measure
  CCV<-function(Nk1,time){ ####Depends on z
    NN<-Nk1+1
    cc<-rep(NA,NN)
    
    cc[1]<- 1/(Nk1)*((0.5*(VCh[1,time]*Tj(0,z[1])+VCh[NN,time]*Tj(0,z[NN])))+
                       sum(VCh[-c(1,NN),time]*Tj(0,z[-c(1,NN)])))
    
    for (i in 2:(Nk1)){
      cc[i]<-2*1/(Nk1)*((0.5*(VCh[1,time]*Tj((i-1),z[1])+VCh[NN,time]*Tj((i-1),z[NN])))+
                          sum(VCh[-c(1,NN),time]*Tj((i-1),z[-c(1,NN)])))
    }
    cc[NN]<-1/(Nk1)*((0.5*(VCh[1,time]*Tj(Nk1,z[1])+VCh[NN,time]*Tj(Nk1,z[NN])))+
                       sum(VCh[-c(1,NN),time]*Tj(Nk1,z[-c(1,NN)])))
    return(cc)
  }
  
  IV<-function(x,Nk1,time){
    # Map x into [-1, 1]
    xx1 <- (2*x - (xop+xned))/(xop-xned)
    
    C<-CCV(Nk1,time) #Creating coef
    len<-Nk1+1
    
    #Evaluate the Chebyshev polynomial
    b2<-rep(0,length(x))
    b1<-rep(0,length(x))
    
    for (i in 1:(len-1)){
      b0<- C[len-(i-1)]+2*xx1*b1-b2
      b2<-b1
      b1<-b0
    }
    f<-C[1]+b1*xx1-b2
    
    return(f)
  }
  
  
  ############################ Parameters
  TT<-(EX+Q)/dt
  end<-dt*TT
  KK<- S(fm,Q,(EX+Q), sigma, lambda)*kp
  
  
  VCh<-matrix(data=NA, nrow=(N+1), ncol=(TT+1))
  CChe<-matrix(data=NA, ncol=(TT+1), nrow=KQ)
  EQChe<- matrix(data=NA, ncol=(TT+1), nrow=KQ)
  
  ######################################## Pre-compute ###################################
  #SWAPTION-Price
  fr<-function(r){
    c<-KK*1
    result<-(1+c)*A(Q+EX,Q+EX+1,sigma,lambda)*exp(-B(Q+EX,Q+EX+1,lambda)*r)
    return((result-1)^2)
  }
  
  rstar<-optimise(fr,c(0.001,0.9))$minimum
  Xstar<-P(Q+EX,Q+EX+1,rstar, lambda, sigma)
  
  Sstar<-100*((1+KK)*ZBC(Q+EX-dt,Q+EX,Q+EX+1,Xstar, sigma, lambda,-0.02))
  Option<-function(x){
    100*((1+KK)*ZBC(Q+EX-dt,Q+EX,Q+EX+1,Xstar, sigma, lambda,x))
  }
  
  
  # Creting the nodal values
  sigmafull<-sqrt(sigma^2/(2*lambda)*(1-exp(-2*lambda*end)))
  meanfull<- alpha(end)+(fm-alpha(0))*exp(-lambda*end)
  xned<- mun(end,0.01)-kk*sigma*sqrt(end)
  xop<- mun(end,0.01)+kk*sigma*sqrt(end)
  z<-cos(pi*seq(0,N)/N)
  xks<-tau(xned,xop,z)
  
  # Calculating expectations of p_j given X=x_k, (k,j)
  sigmap<- sqrt(sigma^2/(2*lambda)*(1-exp(-2*lambda*dt))) #sd of r(t) given r(t-dt)
  mup<-alpha(dt)+(xks-alpha(0))*exp(-lambda*dt) ##Comes from Brigo-Interest rate s. 73
  Glist<-list(Gammat(N,mup,sigmap)) #From period 0 to 1
  times<-seq(dt,(end),dt)
  for (i in 1:(length(times)-1)){ #From period 1 to end
    mup<-alpha(times[i+1])+(xks-alpha(times[i]))*exp(-lambda*dt)
    Glist[[i+1]]<- Gammat(N,mup,sigmap) 
  }
  
  
  
  
  
  ######################### SIMULATION BEGINS ##############################
  startCHE<-Sys.time()
  
  ##### LAST PERIOD
  time<-TT+1
  VCh[,time]<-U(xks,end,end,KK, sigma, lambda)
  EQChe[,time]<-0
  #########
  
  ###### NEXT TO LAST PERIOD
  
  
  time<-time-1
  VCh[,time]<-Option(xks)
  CChe[,time]<-IV(X[,time],N,time)
  
  EQChe[,time]<-CChe[,time]
  
  ab<-Sys.time()
  
  #Before last exercise
  for (j in 1:(exercisestep-2)){
    #print(time-1)
    time<-time-1
    #Coef for VCh(time). From future period.
    LL<-CCV(N,time+1)
    #Computing VCh(time)
    VCh[,time]<-(Glist[[time]]%*%LL)*exp(-dt*xks)
    
    #Computing V(X)
    CChe[,time]<-IV(X[,time],N, time)
    EQChe[,time]<-CChe[,time]
  }
  #Exercise
  time<-time-1
  LL<-CCV(N,time+1)
  
  
  #Computing VCh(time)
  EXhold<-drop((Glist[[time]]%*%LL)*exp(-dt*xks)) #continuation value
  VCh[,time]<-EXhold 
  Cont<-IV(X[,time],N, time)
  print(time)
  VCh[,time]<- pmax(EXhold,U(xks,(time-1)*dt,end,KK, sigma, lambda))
  
  #Computing V(X) and EE
  CChe[,time]<-IV(X[,time],N, time)
  EQChe[,time]<-CChe[,time]
  EXX<- (Cont<=U(X[,time],(time-1)*dt,end,KK, sigma, lambda))
  EQChe[EXX,]<-0
  
  abb<-Sys.time()
  
  #Remaining exercise periods
  for (QW in 1:(EX-1)){
    for (j in 1:(exercisestep-1)){
      
      time<-time-1
      #Coef for VCh(time). From future period.
      LL<-CCV(N,time+1)
      #Computing VCh(time)
      VCh[,time]<-(Glist[[time]]%*%LL)*exp(-dt*xks)
      
      #Computing V(X)
      CChe[,time]<-IV(X[,time],N, time)
      EQChe[,time]<-CChe[,time]
    }
    #Exercise
    time<-time-1
    LL<-CCV(N,time+1)
    
    
    #Computing VCh(time)
    EXhold<-drop((Glist[[time]]%*%LL)*exp(-dt*xks)) #continuation value
    VCh[,time]<-EXhold 
    Cont<-IV(X[,time],N, time)
    print(time)
    VCh[,time]<- pmax(EXhold,U(xks,(time-1)*dt,end,KK, sigma, lambda))
    
    #Computing V(X) and EE
    CChe[,time]<-IV(X[,time],N, time)
    EQChe[,time]<-Cont
    EXX<- (Cont<=U(X[,time],(time-1)*dt,end,KK, sigma, lambda))
    EQChe[EXX,]<-0
  }
  
  #####No exercise
  
  for (b in 1:Q){
    for (j in 1:exercisestep){
      time<-time-1
      #Coef for VCh(time). From future period.
      LL<-CCV(N,time+1)
      #Computing VCh(time)
      a1<-Sys.time()
      VCh[,time]<-(Glist[[time]]%*%LL)*exp(-dt*xks)
      a2<-Sys.time()
      
      #Computing V(X)
      CChe[,time]<-IV(X[,time],N, time)
      a3<-Sys.time()
      EQChe[,time]<-CChe[,time]
    }
  }
  
  
  
  
  ######################################## SIMULERING AF P-measure
  #Parameters
  CCheA<-matrix(data=NA, nrow=KA, ncol=TT) 
  EQCheA<-matrix(data=NA, nrow=KA, ncol=(TT+1))
  
  
  time<-TT+1
  #VCh[,time]<-U(xks,end,end,KK, sigma, lambda)
  EQCheA[,time]<-0
  #########
  
  
  
  time<-time-1
  CCheA[,time]<-IV(XaChe[,time],N,time)
  EQCheA[,time]<-CCheA[,time]
  
  ab<-Sys.time()
  
  #Before last exercise
  for (j in 1:(exercisestep-2)){
    #print(time-1)
    time<-time-1
    
    #Computing V(X)
    CCheA[,time]<-IV(XaChe[,time],N, time)
    EQCheA[,time]<-CCheA[,time]
  }
  #Exercise
  time<-time-1
  LL<-CCV(N,time+1)
  
  
  #Computing VCh(time)
  EXhold<-drop((Glist[[time]]%*%LL)*exp(-dt*xks)) #continuation value
  VCh[,time]<-EXhold 
  Cont<-IV(XaChe[,time],N, time)
  print(time)
  VCh[,time]<- pmax(EXhold,U(xks,(time-1)*dt,end,KK, sigma, lambda))
  
  #Computing V(X) and EE
  CCheA[,time]<-IV(XaChe[,time],N, time)
  EQCheA[,time]<-CCheA[,time]
  EXX<- (Cont<=U(XaChe[,time],(time-1)*dt,end,KK, sigma, lambda))
  EQCheA[EXX,]<-0
  
  #Remaining exercise periods
  for (QW in 1:(EX-1)){
    for (j in 1:(exercisestep-1)){
      
      time<-time-1
      #Computing V(X)
      CCheA[,time]<-IV(XaChe[,time],N, time)
      EQCheA[,time]<-CCheA[,time]
    }
    #Exercise
    time<-time-1
    LL<-CCV(N,time+1)
    
    
    #Computing VCh(time)
    EXhold<-drop((Glist[[time]]%*%LL)*exp(-dt*xks)) #continuation value
    VCh[,time]<-EXhold 
    Cont<-IV(XaChe[,time],N, time)
    print(time)
    VCh[,time]<- pmax(EXhold,U(xks,(time-1)*dt,end,KK, sigma, lambda))
    
    #Computing V(X) and EE
    CCheA[,time]<-IV(XaChe[,time],N, time)
    EQCheA[,time]<-Cont
    EXX<- (Cont<=U(XaChe[,time],(time-1)*dt,end,KK, sigma, lambda))
    EQCheA[EXX,]<-0
    
  }
  
  #####No exercise
  
  for (b in 1:Q){
    for (j in 1:exercisestep){
      time<-time-1
      #Computing V(X)
      CCheA[,time]<-IV(XaChe[,time],N, time)
      EQCheA[,time]<-CCheA[,time]
    }
  }
  
  
  slutCHE<-Sys.time()
  
  
  
  
  ############# Estimation ##################
  print(mean(CChe[,1]))
  #Under P
  PFEChe<-rep(NA,TT+1)
  
  for (i in 1:(TT+1)){
    PFEChe[i]<-quantile(Re(EQCheA[,i]), p=0.99)
  }
  
  mEAChe<- rep(NA,TT+1)
  
  for (i in 1:(TT+1)){
    mEAChe[i]<-mean(EQCheA[,i])
  }
  
  #Under Q
  
  mEQChe<-rep(NA,TT)
  
  for (i in 1:TT){
    mEQChe[i]<-mean(EQChe[,i])
  }
  
  #Discount EQ
  DEQC<-matrix(data=NA, nrow=KQ, ncol=TT+1)
  mDEQC<-rep(NA,TT+1)
  mEQC<-rep(NA,TT+1)
  PFEqC<-rep(NA,TT+1)
  
  AQC<- 0
  i<-1
  DEQC[,i]<- exp(AQC)*EQChe[,i]
  mDEQC[i]<-mean(DEQC[,i])
  PFEqC[i]<- quantile(EQChe[,i],p=0.99)
  mEQC[i]<-mean(EQChe[,i])
  
  
  for (i in 2:(TT+1)){
    AQC<- AQC- 0.5*(X[,i-1]+X[,i])*dt
    DEQC[,i]<- exp(AQC)*EQChe[,i]
    mDEQC[i]<-mean(DEQC[,i])
    mEQC[i]<-mean(EQChe[,i])
    PFEqC[i]<- quantile(EQChe[,i],p=0.99)
  }
  
  #CVA
  dPS<-function(t){
    1-exp(-0.02*t)
  }
  point<-seq(0,end,dt)
  
  
  CVAC<-rep(NA,TT)
  
  
  for (i in 1:(TT)){
    CVAC[i]<- mDEQC[i+1]*(dPS(point[i+1])-dPS(point[i]))
  }
  
  
  
  
  ### Creating output
  
  output<-list()
  
  output[[1]]<-c(mean(CChe[,1]), max(PFEChe), mean(Re(mEAChe[-1])),100*sum(CVAC))
  output[[2]]<-CChe
  output[[3]]<-EQChe
  output[[4]]<-EQCheA
  output[[5]]<-PFEChe
  output[[6]]<-PFEqC
  output[[7]]<-mEAChe
  output[[8]]<- mEQC
  output[[9]]<- slutCHE-startCHE
  output[[10]]<- VCh
  output[[11]]<- xks
  output[[12]]<-Glist
  output[[13]]<- VCh
  
  names(output)<-c("V0, MPFE, EPE, 100CVA", "V", "EQ", "EA", "PFE", "PFEql",  "mEA", "mEQ", "Tid", "VCh", "xks", "Glist", "VCh")
  
  return(output)
  
}
