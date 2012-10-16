as.phylogg <- function(Z,N,tip.label=as.character(1:N)){

#if (is(tip.label)){tip.label=as.character(1:N)}

  edgess<- (matrix(NA,max(Z),2))
  edgel<- (matrix(NA,max(Z),1))

N=as.integer(N)
  
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
 
Zp <- list(edge=edgess, edge.length=edgel,tip.label=tip.label,Nnode=N-1L)
class(Zp) <- "phylo"
return(Zp)
  
}

