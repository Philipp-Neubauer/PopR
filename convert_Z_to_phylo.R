as.phylogg <- function(Z,N,tip.label=as.character(1:N)){

#if (is(tip.label)){tip.label=as.character(1:N)}

Z[Z>N]=Z[Z>N]+1
Z<-rbind(Z,c(N+1,max(Z)+1,1))
  
  edgess<- (matrix(NA,max(Z),2))
  edgel<- (matrix(NA,max(Z),1))

N=as.integer(N+1)
  
z=cbind(Z[,c(2,1)],rep(NA,nrow(Z)))
z[,3]=(N+1):(max(Z)+1)
#z=data.matrix(z)

l = Z[,3]
  
cnt=1;
for (i in 1:nrow(z)){
  if (all(z[i,1:2]<=N))  {
    
    edgess[cnt,] = as.integer(z[i,c(3,1)])
    edgel[cnt,] = l[i]
    edgess[cnt+1,] = as.integer(z[i,c(3,2)])
    edgel[cnt+1,] = l[i]
    cnt=cnt+2
    
  } else if (any(z[i,1:2]<=N))  {
    
    edgess[cnt,] =as.integer( z[i,c(3,1)])
    edgel[cnt,] = l[i]-l[z[,3]==z[i,1]]
    edgess[cnt+1,] = as.integer(z[i,c(3,2)])
    edgel[cnt+1,] = l[i]
    cnt=cnt+2
    
  }else
    { 
     edgess[cnt,] = as.integer(z[i,c(3,1)])
    edgel[cnt,] = l[i]-l[z[,3]==z[i,1]]
     edgess[cnt+1,] = as.integer(z[i,c(3,2)])
    edgel[cnt+1,] = l[i]-l[z[,3]==z[i,2]]
     cnt=cnt+2
    }
}

edgess[edgess[,1]>N,1]= -1L*(edgess[edgess[,1]>N,1]-max(edgess))+N+1L
edgess[edgess[,2]>N,2]= -1L*(edgess[edgess[,2]>N,2]-max(edgess))+N+1L
#edges <- edges[-which(apply(edges,1,function(x){any(is.na(x))})),]
#edges <- edges[order(edges[,1]) ,]
#edgess <-edgess[c(nrow(edgess),1:(nrow(edgess)-1)),]
#edgel <-edgel[c(nrow(edgess),1:(nrow(edgess)-1))]

Zp <- list(edge=edgess, edge.length=edgel,tip.label=c(tip.label,'OUT'),Nnode=N-1L)
class(Zp) <- "phylo"
return(Zp)
  
}
