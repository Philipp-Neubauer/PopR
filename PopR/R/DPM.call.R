DPM.call <-
function(datas=NULL,baseline=NULL,labels=NULL,learn=FALSE,iters=1000,thin=10,np=1,
             typeof='N',path.to.julia=getwd(),call_DPM_path=system.file("exec",package = "PopR"),a.0=NULL,b.0=NULL,ak.0=NULL,bk.0=NULL,v.0=NULL,lambda.0=NULL)             
{
  # initial value for certainty about the mean...
  k.0  = 1
  # initial value for prior mean
  mu.0 = colMeans(datas)

  write.matrix(file='single_priors.csv',c(a.0,b.0,k.0,ak.0,bk.0,v.0,mu.0),sep=',')
  write.matrix(file='matrix_priors.csv',lambda.0,sep=',')
  write.matrix(file='datas.csv',datas,sep=',')

  
  bl=ifelse(learn==F,0,1)
  
  if (learn){
    write.matrix(file='baseline.csv',baseline,sep=',')
    write.matrix(file='labels.csv',labels,sep=',')
  }
  
  if (.Platform$OS.type == "unix")
    {
      exec=file.path(path.to.julia,'./julia')
      command=paste(exec,'-p',1,file.path(call_DPM_path,'call_DPM.jl'),bl,iters,thin,typeof,getwd(),call_DPM_path)
      system(command)
  } else
    {
      exec=file.path(path.to.julia,'/julia/julia.bat')
      command=c('-p',1,file.path(call_DPM_path,'call_DPM.jl'),bl,iters,thin,typeof,getwd(),call_DPM_path)
      system2(exec,command)
  }    
  
 class_ids <- read.csv("source_ids.csv",header=F)
K_record <- read.csv("K_record.csv",header=F)
alpha_record <- read.csv("gammas.csv",header=F)
k_0s <- read.csv("k_0s.csv",header=F)

  list( class_ids= class_ids,K_record=K_record,alpha_record=alpha_record,k_0s=k_0s)
}
