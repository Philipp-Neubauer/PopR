elink.call <- function(class.ids,path.to.julia=getwd(),elink_path=getwd()){

 write.csv(file='class_ids.csv',class.ids)

  if (.Platform$OS.type == "unix")
    {
      exec=file.path(path.to.julia,'./julia')
  } else
    {
      exec=file.path(path.to.julia,'julia/julia.bat')
  }    
 
  command=paste(exec,file.path(elink_path,'elink.jl'),getwd())
  system(command)
  
  link <- read.csv("linkages.csv",header=F)
 # this shouldn't be necessary !
 # link <- link[order(link[,3]),]
  list(tree=link)
}
