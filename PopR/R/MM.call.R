MM.call <-
function(datas=NULL,baseline=NULL,labels=NULL,conditional=FALSE,iters=1000,thin=10,np=1,
             typeof='N',path.to.julia=getwd(),call_MM_path=system.file("exec",package = "PopR"))             
{
  write.csv(file='single_priors.csv',c(a.0,b.0,k.0,ak.0,bk.0,v.0,mu.0))
  write.csv(file='matrix_priors.csv',lambda.0)
  write.csv(file='datas.csv',datas)

  
  cond=ifelse(conditional==F,0,1)
  
  
    write.csv(file='baseline.csv',baseline)
    write.csv(file='labels.csv',labels)
  
  
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
