###### Functions regarding the pricing of Hull-White swaption and further assumption made in the Qian-paper

#Swap rate
S<-function(x,n,N, sigma, lambda){
  a<-0
  for (i in n:N){
    a<- a+P(n,i+1,x, sigma, lambda)
  }
  return((1-P(n,N+1,x, sigma, lambda))/a)
}


#Payoff if option is exercised
U<-function(x,n,N, K, sigma, lambda){
  N0<-100
  a<-0
  for (i in n:N){
    a<-a+P(n,i+1,x, sigma, lambda)
  }
  
  Ut<- N0*a*pmax(-S(x,n,N, sigma, lambda)+K,0)
  return(Ut)
}


#ZCB from the market
Pm<-function(TT){
  exp(-0.01*TT)
}


### ZCB
B<-function(t,TT, lambda){
  1/lambda* (1-exp(-lambda*(TT-t)))
}

A<-function(t,TT, sigma, lambda){
  Pm(TT)/Pm(t)*exp(B(t,TT, lambda)*fm-sigma^2/(4*lambda)*(1-exp(-2*lambda*t))*B(t,TT,lambda)^2)  
}

P<-function(t,TT,x, sigma, lambda){
  A(t,TT, sigma, lambda)*exp(-B(t,TT, lambda)*x) #(Mistake in paper (QIAN))
}

ZBC<-function(t,TT,S,X,sigma,lambda,x){ #Bringo 'smile volatility'
  sigmap<-sigma*sqrt((1-exp(-2*lambda*(TT-t)))/(2*lambda))*B(TT,S,lambda)
  h<- 1/sigmap *log(P(t,S,x,sigma,lambda)/(P(t,TT,x,sigma,lambda)*X))+sigmap/2
  return(P(t,S,x,sigma,lambda)*pnorm(h)-X*P(t,TT,x, sigma, lambda)*pnorm(h-sigmap))
}


#Simulating of the rate
alpha<-function(t){ 
  0.01+sigma^2/(2*lambda^2)*(1-exp(-lambda*t))^2 #Comes from Brigo-Interest rate s. 73
}

r<-function(t,TT,x,k){
  m  <- alpha(TT)+(x-alpha(t))*exp(-lambda*(TT-t)) ##Comes from Brigo-Interest rate s. 73
  s2 <- sigma^2*(1-exp(-2*lambda*(TT-t)))/(2*lambda) #Comes from Brigo-Interest rate s. 73
  Z <- rnorm(k,m,sqrt(s2))
  return(Z)
}
