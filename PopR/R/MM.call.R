MM.call <-
function(datas=NULL,baseline=NULL,labels=NULL,conditional=FALSE,iters=1000,thin=10,np=1,
             typeof='N',path.to.julia=getwd(),call_MM_path=system.file("exec",package = "PopR"),v.0=NULL,lambda.0=NULL)             
{
  # initial value for certainty about the mean...legacy
  k.0  = 1
  a.0=0
  b.0=0
  ak.0=0
  bk.0=0
  # initial value for prior mean
  mu.0 = colMeans(datas)
  write.matrix(file='single_priors.csv',c(a.0,b.0,k.0,ak.0,bk.0,v.0,mu.0),sep=',')
  write.matrix(file='matrix_priors.csv',lambda.0,sep=',')
  write.matrix(file='datas.csv',datas,sep=',')

  
  cond=ifelse(conditional==F,0,1)
  
  
   write.matrix(file='baseline.csv',baseline,sep=',')
    write.matrix(file='labels.csv',labels,sep=',')
  
  
  if (.Platform$OS.type == "unix")
    {
      exec=file.path(path.to.julia,'./julia')
      command=paste(exec,'-p',1,file.path(call_MM_path,'call_MM.jl'),cond,iters,thin,typeof,getwd(),call_MM_path)
      system(command)
  } else
    {
      exec=file.path(path.to.julia,'/julia/julia.bat')
      command=c('-p',1,file.path(call_MM_path,'call_MM.jl'),cond,iters,thin,typeof,getwd(),call_MM_path)
      system2(exec,command)
  }    
  
 class_ids <- read.csv("source_ids.csv",header=F)
proportions <- read.csv("proportions.csv",header=F)
post_probas <- read.csv("post_probas.csv",header=F)


  list( class_ids= class_ids,proportions=proportions,post_probas=post_probas)
}
