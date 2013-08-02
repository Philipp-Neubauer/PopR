DPM.call <-
function(datas=NULL,baseline=NULL,labels=NULL,learn=FALSE,iters=1000,thin=10,np=1,
             typeof='N',path.to.julia=getwd(),call_DPM_path=system.file("exec",package = "PopR"))             
{
  write.csv(file='single_priors.csv',c(a.0,b.0,k.0,ak.0,bk.0,v.0,mu.0))
  write.csv(file='matrix_priors.csv',lambda.0)
  write.csv(file='datas.csv',datas)

  
  bl=ifelse(learn==F,0,1)
  
  if (learn){
    write.csv(file='baseline.csv',baseline)
    write.csv(file='labels.csv',labels)
  }
  
  if (.Platform$OS.type == "unix")
    {
      exec=file.path(path.to.julia,'./julia')
      command=paste(exec,'-p',np,file.path(call_DPM_path,'call_DPM.jl'),bl,iters,thin,typeof,getwd(),call_DPM_path)
      system(command)
  } else
    {
      exec=file.path(path.to.julia,'/julia/julia.bat')
      command=c('-p',np,file.path(call_DPM_path,'call_DPM.jl'),bl,iters,thin,typeof,getwd(),call_DPM_path)
      system2(exec,command)
  }    
  
 class_ids <- read.csv("source_ids.csv",header=F)
K_record <- read.csv("K_record.csv",header=F)
alpha_record <- read.csv("gammas.csv",header=F)
k_0s <- read.csv("k_0s.csv",header=F)

  list( class_ids= class_ids,K_record=K_record,alpha_record=alpha_record,k_0s=k_0s)
}
