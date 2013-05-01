########### Toy data example ############

# DPM toolbox by P.Neubauer -  distributed under GNU GPL licence V3.

require(MASS)
require(ape)
source("Julia_call_function.R")
source("elink.call.R")
source("convert_Z_to_phylo.R")
source("prior_match.R")
setwd("/home/philbobsqp/Work//Projects/Papersampler")
## task 1 - estiamte correct number of sources without a baseline

#  multi-normal source distributions - in an ideal world...

num.sources = 5  # 'true' number of sources
num.elements = 8 # number of elements
num.per.source = 20 # individuals per source

sep = 10 # separation of means

means = mvrnorm(num.sources,rep(0,num.elements),diag(rep(sep,num.elements)))

datas=matrix(NA,num.per.source*num.sources,num.elements)
label=rep(NA,num.per.source*num.sources)
a=0
for (i in 1:num.sources){
    datas[(a*num.per.source+1):(a*num.per.source+num.per.source),] <-mvrnorm(num.per.source,means[i,],diag(rep(1,num.elements)))
    label[(a*num.per.source+1):(a*num.per.source+num.per.source)] <- i
    a=a+1
}
label=sort(rep(1:num.sources,num.per.source))
scores = princomp(datas,label)$scores
#pdf('PCAplot.pdf',colormodel='cmyk')
plot(scores[,1],scores[,2],col=label,xlab='LD1',ylab='LD2')
#dev.off()
##########################
### now use the DPM ######
##########################

data.DPM = t(apply(datas,1,function(x){x-colMeans(datas)}))

# prior for gamma - set to uninformative following Dorazio 2009
n=num.per.source*num.sources

# defaults to uniform prior, can take a vector or probability mass, or arguments poisson, negbin,norm,lnorm, tha last 4 have parameters mu and var/rate
ab = get_prior_ab(n,"poisson",5)


a.0 = ab[[1]]
b.0 = ab[[2]]
a.0 
b.0

#Inv-Wishart prior

# degrees of freedom - at least num.elements+1
# higher v.0 equate to a more informative prior
v.0  = num.elements+1

# prior covariance - try a number of options - this shows how sensitive the
# model/data combination is to prior choice...

# play with the variance to see how the variance and source separation
# interact (use same var for each element.) - for good source separation
# this should not be very important i.e. the prior should not influence the
# number of sources. This changes when sources are not easily identifyable

vars = 1 # 
#adjust prior by degrees of freedom to get the right expected value
lambda.0 = diag(rep(vars*(v.0-num.elements),num.elements))

# prior mean for k_0 is ak.0*bk.0
ak.0 = 1
bk.0 = 1

ak.0*bk.0

# initial value for certainty about the mean...
k.0  = 1
# initial value for prior mean
mu.0 = colMeans(data.DPM)

# number of iterations per processor
num.iters=1000
# thinning of the marcov chain
thin=1

# number of parallel processes (chains) to run (recommended to keep at <= core of CPU) set to np + 1 since one instance only calls and summarizes the work...
np=4+1

# total number of kept iterations
niter=(np-1)*num.iters/thin
burnin = 100  # number of (kept!) iterations to discard


############## Run the sampler ##############

# there will most likely be no output on the terminal in windows until the very end. 
#this works better in Linux where progress is displayed continously - K^+ is the estimated number of sources

output = DPM.call(datas=data.DPM,iters=num.iters,thin=thin,np=np, path.to.julia='/home/philbobsqp/Work/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)


############## Analyse output ###############

# check if the markov chain for number of sources has converged and is mixing:
plot(classes[burnin:niter,1]) # number of sources
plot(output$alpha_record[burnin:niter,1]) #concentration parameter 
plot(output$k_0[burnin:niter,1]) # prior k_0

# display the histogram of the number of sources

bins = (min(classes[burnin:niter,1])-0.5):(max(classes[burnin:niter,1])+0.5) # histogram bins

hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)

# now create the exact linkage tree and display

S=class.id[,burnin:niter]
Z = elink.call(S, path.to.julia='/home/philbobsqp/Work/julia')$tree
Zp <- as.phylogg(Z,n,rep('o',n))

#pdf('./Plots/tree very easy example.pdf')
plot.phylo(reorder(Zp, order = "c"),edge.width=2,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(label,0),type='f')
#dev.off()


hc = hclust(dist(data.DPM))
hcp = as.phylo(hc)
hcp$tip.label=rep('o',N)

pdf('./Plots/clustering example.pdf')
plot.phylo(reorder(hcp, order = "c"),edge.width=2,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(label,0),type='f')
dev.off()

pdf('hardest5 example.pdf',width=12, height=3,colormodel='cmyk')
par(mfrow=c(1,4))
par(mar=c(5,4,2,1)+0.1)
plot(scores[,1],scores[,2],col=label,xlab='LD1',ylab='LD2')
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(hcp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length))),tip.color=c(label),type='f')
par(mar=c(5,4,2,0)+0.1)
hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(Zp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(label,0),type='f')
dev.off()

#### playing with the covariance shows the sensitivity to the prior
#### component variance one needs to define a source in terms of variance
#### about a mean to be able to estimate its number. I would personally not
#### recommend trying to estimate the number of sources without any prior information as it must fail (...the
#### sum of two normals is a normal...or inversely, the normal distribution is infinitiely divisible...)

## 
###################################
### Using the DPM with a baseline #
#########################################

# same as before, set up the data
num.sources = 4  # sampled number of sources,
extra.sources = 2 # extra baseline
num.elements = 5 # number of elements
num.per.source = rep(30,num.sources) # number of samples per source
num.per.source = c(num.per.source,sample.int(30,extra.sources)+5) # number of samples per extra.sources
as=num.sources+extra.sources

sep = 8 # separation of means
means = mvrnorm(as,rep(0,num.elements),diag(rep(sep,num.elements)))

data=matrix(NA,sum(num.per.source),num.elements)
label=rep(NA,sum(num.per.source))
a=0
for (i in 1:as){
    a=a+1
    data[a:(a+num.per.source[i]-1),] <-mvrnorm(num.per.source[i],means[i,],diag(rep(1,num.elements)))
    label[a:(a+num.per.source[i]-1)] <- i
    a=a+num.per.source[i]-1
}


scores = princomp(data)$scores

pdf('../easy example fix.pdf',colormodel='cmyk')
plot(scores[,1],scores[,2],t='n',xlab='PCA1',ylab='PCA2')

points(scores[1:sum(num.per.source[1:num.sources]),1],scores[1:sum(num.per.source[1:num.sources]),2],col=label[1:sum(num.per.source[1:num.sources])],pch=21,bg=label[1:sum(num.per.source[1:num.sources])])

points(scores[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),1],scores[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),2],col=label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)],pch=23,bg=label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)])

dev.off()


##########################
### now use the DPM ######
##########################

baseline = data[1:sum(num.per.source[1:num.sources]),]

# add 1/4th of the samples from the baseline to the mixed set

ixs = sample.int(sum(num.per.source[1:num.sources]),round(sum(num.per.source[1:num.sources])/4),1)

mixed = rbind(baseline[ixs,],data[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),])
baseline=baseline[-ixs,]

bix = 1:sum(num.per.source[1:num.sources])
bix=bix[-ixs]

mixedlabels=label[c(ixs,(sum(num.per.source[1:num.sources])+1):sum(num.per.source))]
baselabels = label[bix]

mixed = t(apply(mixed,1,function(x){x-colMeans(baseline)}))
baseline = t(apply(baseline,1,function(x){x-colMeans(baseline)}))


# prior for gamma - set to uninformative following Dorazio 2009
n=sum(num.per.source[(num.sources+1):(num.sources+extra.sources)])
ab = get_prior_ab(n,'norm',14,1)

a.0 = ab[[1]]
b.0 = ab[[2]]
a.0
b.0

#Inv-Wishart prior

# degrees of freedom - at least num.elements+1
# higher v.0 equates to a more informative prior
v.0  = num.elements+1

# prior covariance - try a number of options - this shows how sensitive the
# model/data combination is to prior choice...

# play with the variance to see how the variance and source separation
# interact (use same var for each element.) - for good source separation
# this should not be very important i.e. the prior should not influence the
# number of sources. This changes when sources are not easily identifyable

vars = by(baseline,baselabels,cov) # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
var=0
for (i in 1:num.sources){
  var=var+(1/num.sources)*solve(vars[[i]])
}
lambda.0 = solve(var)
#lambda.0=diag(0.01,num.elements)

# prior for k_0 - small values are uninformative, but may lead to very poor mixing and numerical instability.
ak.0 = 1
bk.0 = 1
# initial k_0...
k.0  = 1
# prior mean
mu.0 = colMeans(baseline)

# number of iterations per processor
num.iters=2000

# numebr of parallel processing jobs - set to a max of number of cores of the CPU+1(the +1 just calls the cumputing instances)
np=2+1
# thinning of the marcov chain
thin=2
burnin = 100  # number of (kept!) iterations to discard

# total number of kept iterations
niter=(np-1)*num.iters/thin

# if julia is not installed and in the path, or is moved to a different directory, the path to the executeable must be provided, else the working directory is taken as default
output =DPM.call(datas=mixed,learn=T,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels,path.to.julia='/home/philbobsqp/Work/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)

# check if the markov chain for number of sources has converged and is mixing:
plot(classes[burnin:niter,1]) # number of sources
plot(output$alpha_record[burnin:niter,1]) #concentration parameter 
plot(output$k_0[burnin:niter,1]) # prior k_0

# display the histogram of the number of sources

burnin = 100  # number of (kept!) iterations to discard
bins = (min(classes[burnin:niter,1])-0.5):(max(classes[burnin:niter,1])+0.5) # histogram bins

hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)

# now create the exact linkage tree and display - for this one the first
# num.sources leafs are the baseline
S.fix=class.id[,burnin:niter]
for (i in num.sources:1){
S.fix= rbind(rep(i,niter-burnin),S.fix)}

Z2 = elink.call(S.fix,'/home/philbobsqp/Work/julia')$tree
N=length(mixedlabels)
Zp <- as.phylogg(Z2,N+num.sources,c(rep('X',num.sources),rep('o',N)))


plot.phylo(reorder(Zp, order = "c"),edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(1:num.sources,mixedlabels,0),type='f')



hc = hclust(dist(data.DPM))
hcp = as.phylo(hc)
hcp$tip.label=rep('o',N)


pdf('fix simulation.pdf',width=5, height=5,colormodel='cmyk')
par(mar=c(5,4,2,1)+0.1)
plot(scores[,1],scores[,2],t='n',xlab='PCA1',ylab='PCA2')
points(scores[1:sum(num.per.source[1:num.sources]),1],scores[1:sum(num.per.source[1:num.sources]),2],col=label[1:sum(num.per.source[1:num.sources])],pch=21,bg=label[1:sum(num.per.source[1:num.sources])])
points(scores[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),1],scores[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),2],col=label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)],pch=23,bg=label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)])
dev.off()

pdf('fix low var example.pdf',width=6, height=3,colormodel='cmyk')
par(mfrow=c(1,2))
par(mar=c(5,4,2,0)+0.1)
hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(Zp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(1:num.sources,mixedlabels,0),type='f')
dev.off()

# lda classification succes
baselda <- lda(baseline,baselabels)
mean(apply(predict(baselda,mixed)$posterior,1,which.max)==mixedlabels)

# compare against mixture models

# conditional analysis
cond.output =MM.call(datas=mixed,conditional=T,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels,path.to.julia='/home/philbobsqp/Work/julia')

cond.class.id = as.data.frame(cond.output$class_id)
cc.fix=cond.class.id[,burnin:niter]
for (i in num.sources:1){
cc.fix= rbind(rep(i,niter-burnin),cc.fix)}

Z = elink.call(cc.fix,path.to.julia='/home/philbobsqp/Work/julia')$tree
N=length(mixedlabels)
Zp <- as.phylogg(Z,N+num.sources,c(rep('X',num.sources),rep('o',N)))


plot.phylo(reorder(Zp, order = "c")
           ,edge.color=c(rep(1,length(Zp$edge.length)-1),0)
           ,tip.color=c(1:num.sources,mixedlabels,0)
           ,type='f')


# unconditional analysis
uncond.output =MM.call(datas=mixed,conditional=F,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels,path.to.julia='/home/philbobsqp/Work/julia')

ucond.class.id = as.data.frame(uncond.output$class_id)
uc.fix=ucond.class.id[,burnin:niter]
for (i in num.sources:1){
uc.fix= rbind(rep(i,niter-burnin),uc.fix)}

Z = elink.call(uc.fix,path.to.julia='/home/philbobsqp/Work/julia')$tree
N=length(mixedlabels)
Zp <- as.phylogg(Z,N+num.sources,c(rep('X',num.sources),rep('o',N)))


plot.phylo(reorder(Zp, order = "c")
           ,edge.color=c(rep(1,length(Zp$edge.length)-1),0)
           ,tip.color=c(1:num.sources,mixedlabels,0)
           ,type='f')
