get_Pc <- function(inds1,inds2,phylo){
  
  return(1-mean(cophenetic.phylo(phylo)[inds1,inds2])/2)
  
}