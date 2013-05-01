get_prior_ab <- function(n,g='uniform',mu=NULL,var=NULL){
  if (!is.numeric(g)){
  if(g=='uniform'){
    g <- rep(1/n,n)
  } else if (g=='poisson') {
    g <- dpois(1:n,mu)
  } else if (g=='negbin') {
    g <- dnbinom(1:n,mu,var)
  } else if (g=='lnorm') {
    g <- dlnorm(1:n,mu,var)
  } else if (g=='norm'){
    g <- dnorm(1:n,mu,var)
  } else {stop('Distribution need to be one of the following: uniform (default),poisson,negbin,lnorm,norm')}
}
g=g/sum(g)

s_recursive <- function(n,k,s_rec) {
    if(k==1 & n==1){
        return(1)}else if(k==1){return(  -(n-1)*s_rec[n-1,k])
        }else{ 
      return(s_rec[n-1,k-1] - (n-1)*s_rec[n-1,k])
    }
  }

  # get stirling numbers
  snk = matrix(,n,n);
  for (i in 1:n){
    for (j in 1:(i)){
      if (i==j) snk[i,j] =1 else snk[i,j] = (s_recursive(i,j,snk) )
    }
}
snk=abs(snk)
  
  # KL distance
  KL_dist <- function(par,others){
    
    a <-par[1]
    b <-par[2]
    n <- others$n
    g <- others$g
    snk <- others$snk
    
    #integral
    pi_fun <- function(alpha,otros){
      a <-otros$a
      b <-otros$b
      n <- otros$n
      k <- otros$k
      exp((k+a-1)*log(alpha)-b*alpha+lgamma(alpha)-lgamma(alpha+n))
      
    }
    
    KL = 0;pi_k=vector(,n)
    for (k in 1:n){
      nint <- integrate(pi_fun,lower=0,upper=Inf,otros=list(a=a,b=b,k=k,n=n),stop.on.error =F,abs.tol=0.)$value
      pi_k[k] <- ((b^a*snk[n,k])/gamma(a))*nint
      
    }
    #pi_k=pi_k/sum(pi_k)
    KL <- sum(g*log(g/pi_k))
   # cat(KL,'\n')
  
    return(KL)
  }
  
  res <- optim(par=c(0.5,0.5),KL_dist,others=list(n=n,g=g,snk=snk),lower=c(1E-3,1E-4),upper=c(100,100),method='L-BFGS-B')  
  
  a <- res$par[1]
  b <- res$par[2]
  
#   # same as within the KL function, just for plotting
    
    #integral
    pi_fun <- function(alpha,otros){
      a <-otros$a
      b <-otros$b
      n <- otros$n
      k <- otros$k
      exp((k+a-1)*log(alpha)-b*alpha+lgamma(alpha)-lgamma(alpha+n))
      
    }
    
    pi_k=vector(,n)
    for (k in 1:n){
      nint <- integrate(pi_fun,lower=0,upper=Inf,otros=list(a=a,b=b,k=k,n=n),stop.on.error =T,abs.tol=0.)$value
      pi_k[k] <- (((b^a)*snk[n,k])/gamma(a))*nint
      
    }
  
     pi_k=pi_k/sum(pi_k)
     plot(1:n,pi_k,pch=16,col=4,xlab='number of sources',ylab='prior probability mass')
     points(1:n,g,pch=17,col=3)
  legend(n-n/3,(max(c(pi_k,g))),pch=16:17,legend=c("induced prior","original prior"),col=4:3)
  list(a=a,b=b)
  
}
