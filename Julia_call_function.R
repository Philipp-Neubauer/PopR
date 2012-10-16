# bl = baseline present ?
# iters = number of MCMC iterations
# thin = thinning of the markov chain(s)
# np = number of processes = number of chains
# !!! np must be less or equal to number of available proceesors/cores (e.g., 2 on a dual core machine, or 3 on a quad core machine)
# !!! WARNING setting np=number of cores  may take away all or most of your computing resources

DPM.call <- function(baseline=F,iters=10000,thin=10,np=1,
             path.to.julia=getwd(),call_DPM_path=getwd())             
{

  write.csv(file='single_priors.csv',c(a.0,b.0,k.0,v.0,mu.0))
  write.csv(file='matrix_priors.csv',lambda.0)
  write.csv(file='datas.csv',data.DPM)

  if (baseline){
    write.csv(file='baseline.csv',baseline)
    write.csv(file='labels.csv',labels)
  }
  
  if (.Platform$OS.type == "unix")
    {
      exec=file.path(path.to.julia,'./julia')
  } else
    {
      exec=file.path(path.to.julia,'julia')
  }    

bl=ifelse(baseline==F,0,1)
  
  command=paste(exec,'-p',np,file.path(call_DPM_path,'call_DPM.jl'),bl,iters,thin,getwd())
  system(command)
  
 class_ids <- read.csv("source_ids.csv",header=F)
K_record <- read.csv("K_record.csv",header=F)
alpha_record <- read.csv("gammas.csv",header=F)
k_0s <- read.csv("k_0s.csv",header=F)

  list( class_ids= class_ids,K_record=K_record,alpha_record=alpha_record,k_0s=k_0s)
}



