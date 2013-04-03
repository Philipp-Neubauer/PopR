get_prior_ab <- function(n,g='uniform',pars=NULL){
  
  if(g=='uniform'){
    g <- rep(1/n,n)
  } else if (g=='poisson') {
    g <- dpois(1:n,pars)
  } else if (g=='negbin') {
    g <- dnbinom(1:n,pars)
  }
  
#   s_nk <- function(n,k){
#     
#     kt=0;
#     for(m in 0:(n-k))  {
#       rt = 0
#       for (j in 0:m) rt=rt+(((-1)^j*j^(n-k+m))/(factorial(j)*factorial(m-j)))
#       kt <- kt + (1/((n+m)*factorial(n-k-m)*factorial(n-k+m)))*rt
#     }
#     return((factorial(2*n-k)/factorial(k-1))*kt)
#   }
#   
#   snk <- matrix(,n,n)
#   for (i in 1:n){
#     for (j in 1:i){
#       snk[i,j]<-abs(s_nk(i,j))
#     }
#   }
  
  KL_dist <- function(par,others){
    # integral 
    a <-par[1]
    b <-par[2]
    n <- others$n
    g <- others$g
    #snk <- others$snk
    
    pi_fun <- function(alpha,otros){
      a <-otros$a
      b <-otros$b
      n <- otros$n
      k <- otros$k
      exp((k+a-1)*log(alpha)-b*alpha+lgamma(alpha)-lgamma(alpha+n))
      
    }
    
    KL = 0
    for (k in 1:n){
      nint <- integrate(pi_fun,lower=0,upper=Inf,otros=list(a=a,b=b,k=k,n=n),stop.on.error =F,abs.tol=0.)$value
      #nint <- quadgr(pi_fun,0,Inf,otros=list(a=a,b=b,k=k,n=n),tol=0.)$value
      pi_k <- ((b^a)/gamma(a))*nint
      KL <- KL + g[k]*log(g[k]/pi_k)
      #KL <- KL +log(pi_k)
      #cat(log(pi_k),'\n')
    }
    #KL <- -log(n)-KL/n 
    return(KL)
  }
  
  res <- optim(par=c(0.47,0.006),KL_dist,others=list(n=n,g=g),lower=c(1E-3,1E-4),upper=c(100,100),method='L-BFGS-B')  
  
  a <- res$par[1]
  b <- res$par[2]
  
#   pi_fun <- function(alpha,a,b,k,n){exp((k+a-1)*log(alpha)-b*alpha+lgamma(alpha)-lgamma(alpha+n))}
#   
#   pi_k <- vector(,n)
#   for (k in 1:n){
#     nint <- integrate(pi_fun,lower=0,upper=Inf,a=a,b=b,k=k,n=n,stop.on.error =F, abs.tol = 0.)
#     pi_k[k] <- (((b^a))/gamma(a))*nint$value
#   }
#   
#   plot(1:n,pi_k,xlab='number of sources',ylab='prior probability mass')
#   
  list(a=a,b=b)
  
}