########### Toy data example ############

# DPM toolbox by P.Neubauer -  distributed under GNU GPL licence V3.

require(MASS)
require(ape)
source("Julia_call_function.R")
source("elink.call.R")
source("convert_Z_to_phylo.R")
source("prior_match.R")
setwd("/home/philbobsqp/Work//Projects/Papersampler")

# set color palette to colorblind friendly colors
cbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#FF796B","#999999", "#F0E442", "#0072B2", "#D55E00")

## task 1 - estiamte correct number of sources without a baseline
## USE WITH CAUTION, THIS IS A DIFFICULT< IF NOT IMPOSSIBLE TASK IN MANY REALISTIC SETTINGS 
## WITHOUT ANY PRIOR INFO

#  multi-normal source distributions - in an ideal world...

num.sources = 4  # 'true' number of sources
num.elements = 5 # number of elements
num.per.source = 30 # individuals per source

sep = 24 # separation of means

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
scores = princomp(datas)$scores
#pdf('PCAplot.pdf',colormodel='cmyk')
plot(scores[,1],scores[,2],col=cbPalette[label],xlab='PC1',ylab='PC2',pch=16)
#dev.off()
##########################
### now use the DPM ######
##########################

# use original data
data.DPM = t(apply(datas,1,function(x){x-colMeans(datas)}))
num.elements=ncol(data.DPM)

# OR use first num.elements principal components to overcome dimensionality problems: 
# NOTE, this may throw away useful info --- 
# summary(princomp(datas))
# num.elements=2
# data.DPM = t(apply(scores[,1:num.elements],1,function(x){x-colMeans(scores[,1:num.elements])}))


# prior for gamma - set to matching prior following Dorazio 2009
n=num.per.source*num.sources

# defaults to uniform prior, can take a vector or probability mass, or arguments poisson, negbin,norm,lnorm, tha last 4 have parameters mu and var/rate
# best case scenario, we have a vargue idea that the number of sources in the data is 5
ab = get_prior_ab(n,"poisson",4) # for uniform omit last two arguments to the function


a.0 = ab[[1]]
b.0 = ab[[2]]
a.0 
b.0

#Inv-Wishart prior

# degrees of freedom - at least num.elements+1
# higher v.0 equate to a more informative prior for source covariance
v.0  = num.elements+20

vars = 1 # best case scenario - we actually know the (co-)variance...this is not the case in practice

#adjust cov prior by degrees of freedom to get the right expected value
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
np=2+1

# total number of kept iterations
niter=(np-1)*num.iters/thin
burnin = 100  # number of (kept!) iterations to discard


############## Run the sampler ##############

# there will most likely be no output on the terminal in windows until the very end. 
#this works better in Linux (OSX?)where progress is displayed continously - K^+ is the estimated number of sources

output = DPM.call(datas=data.DPM,iters=num.iters,thin=thin,np=np, path.to.julia='/home/philbobsqp/Work/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)


############## Analyse output ###############

# check if the markov chain for number of sources has converged and is mixing:
lpc=niter/(np-1)
keeps = rep((burnin+1):(niter/(np-1)),np-1)+rep(0:(np-2)*lpc,each=niter/2-burnin)

plot(classes[keeps,1]) # number of sources
plot(output$alpha_record[keeps,1]) #concentration parameter 
plot(output$k_0[keeps,1]) # prior k_0

# display the histogram of the number of sources

bins = (min(classes[keeps,1])-0.5):(max(classes[keeps,1])+0.5) # histogram bins

hist(classes[keeps,1],bins,col='grey',xlab='number of sources',main='',freq=F)

# now create the exact linkage tree and display

S=class.id[,keeps]
Z = elink.call(S, path.to.julia='/home/philbobsqp/Work/julia')$tree
Zp <- as.phylogg(Z,n,rep('o',n))

#pdf('./Plots/tree very easy example.pdf')
plot.phylo(reorder(Zp, order = "c"),edge.width=2,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(cbPalette[label],0),type='f')
#dev.off()

hc = hclust(dist(data.DPM))
hcp = as.phylo(hc)
hcp$tip.label=rep('o',n)

pdf('hardest PCAPD Poisprior example.pdf',width=12, height=3)
par(mfrow=c(1,4))
par(mar=c(5,4,2,1)+0.1)
plot(scores[,1],scores[,2],col=cbPalette[label*1.8],xlab='PCA1',ylab='PCA2')
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(hcp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(hcp$edge.length))),tip.color=c(cbPalette[label*1.8]),type='f')
par(mar=c(5,4,2,0)+0.1)
hist(classes[keeps,1],bins,col='grey',xlab='number of sources',main='',freq=F)
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(Zp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(cbPalette[label*1.8],0),type='f')
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
num.per.source = rep(50,num.sources) # number of samples per source
num.per.source = c(num.per.source,sample.int(10,extra.sources)+5) # number of samples per extra.sources
as=num.sources+extra.sources

sep = 4 # separation of means
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

#linear discriminant analysis for plotting
scores = lda(data[1:sum(num.per.source[1:num.sources]),],label[1:sum(num.per.source[1:num.sources])])$scaling

#pdf('../easy example fix.pdf',colormodel='cmyk')
plot(data%*%scores[,1],data%*%scores[,2],t='n',xlab='PCA1',ylab='PCA2')

# poitns are baseline sources
plot(data[1:sum(num.per.source[1:num.sources]),]%*%scores[,1],data[1:sum(num.per.source[1:num.sources]),]%*%scores[,2],col=cbPalette[label[1:sum(num.per.source[1:num.sources])]],pch=21,bg=cbPalette[label[1:sum(num.per.source[1:num.sources])]],xlab='LDA1',ylab='LDA2')
# diamonds are extra-baseline
points(data[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),]%*%scores[,1],data[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),]%*%scores[,2],col=cbPalette[label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)]],pch=23,bg=cbPalette[label[(sum(num.per.source[1:num.sources])+1):sum(num.per.source)]])

#dev.off()


##########################
### now use the DPM ######
##########################

baseline = data[1:sum(num.per.source[1:num.sources]),]

# add 1/4th of the samples from the baseline to the mixed set

ixs = sample.int(sum(num.per.source[1:num.sources]),round(sum(num.per.source[1:num.sources])/3),1)

mixed = rbind(baseline[ixs,],data[(sum(num.per.source[1:num.sources])+1):sum(num.per.source),])
baseline=baseline[-ixs,]

bix = 1:sum(num.per.source[1:num.sources])
bix=bix[-ixs]

mixedlabels=label[c(ixs,(sum(num.per.source[1:num.sources])+1):sum(num.per.source))]
baselabels = label[bix]

mixed = t(apply(mixed,1,function(x){x-colMeans(baseline)}))
baseline = t(apply(baseline,1,function(x){x-colMeans(baseline)}))

## backups for testing

#mix_o = mixed 
#base_0 = baseline

mixed_t = mixed%*%scores
baseline_t = baseline%*%scores
#num.elements=ncol(mixed)

# rounds are baseline, triangles are mixed sample

plot(baseline_t[,1],baseline_t[,2],col=cbPalette[baselabels],pch=16)
points(mixed_t[,1],mixed_t[,2],col=cbPalette[mixedlabels],pch=17)

#mixed = mix_o
#baseline=base_0 
num.elements=ncol(mixed)
## done


# prior for gamma - set to something more informative now
n = length(mixedlabels)
ab = get_prior_ab(n,'poisson',4)

a.0 = ab[[1]]
b.0 = ab[[2]]
a.0
b.0

#Inv-Wishart prior

# degrees of freedom - at least num.elements+1
# higher v.0 equates to a more informative prior
v.0  = num.elements+1

# prior covariance 
vars = by(baseline,baselabels,cov) 
var=0
for (i in 1:num.sources){
  var=var+(1/num.sources)*solve(vars[[i]])
}
lambda.0 = solve(var)

# prior for k_0 - small values are uninformative, but may lead to very poor mixing and numerical instability.
ak.0 = 1
bk.0 = 1
# initial k_0...is estimated so no need to change
k.0  = 1
# initial prior mean...is estimated so no need to change
mu.0 = colMeans(baseline)

# numebr of parallel processing jobs - set to a max of number of cores of the CPU+1(the +1 just calls the cumputing instances)
np=1+1

burnin = 100  # number of (kept!) iterations to discard

# number of iterations per processor
num.iters=1000+burnin
# thinning of the Markov chain
thin=2
# total number of kept iterations
niter=(np-1)*num.iters/thin

# if julia is not installed and in the path, 
# or is moved to a different directory, the 
#path to the executeable must be provided (and changed in the line below), 
#else the working directory is taken as default
outputs =DPM.call(datas=mixed,learn=T,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels,path.to.julia='/home/philbobsqp/Work/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(outputs$class_id)
# and just the number of sources per iteration
classes = as.data.frame(outputs$K_record)

# check if the markov chain for number of sources has converged and is mixing:
lpc=niter/(np-1)
keeps = rep((burnin+1):(niter/(np-1)),np-1)+rep(0:(np-2)*lpc,each=niter/2-burnin)

plot(classes[keeps,1]) # number of sources
plot(outputs$alpha_record[keeps,1]) #concentration parameter 
plot(outputs$k_0[keeps,1]) # prior k_0

# display the histogram of the number of sources

bins = (min(classes[keeps,1])-0.5):(max(classes[keeps,1])+0.5) # histogram bins

hist(classes[keeps,1],bins,col='grey',xlab='number of sources',main='',freq=F)

# now create the exact linkage tree and display - for this one the first
# num.sources leafs are the baseline - S.fix just appends the 
#(collapsed - the co-assignment is allways 1) baseline to the sampled source ids in class.id
S.fix=class.id[,keeps]
for (i in num.sources:1){
S.fix= rbind(rep(i,niter-burnin),S.fix)}

# make tree
Z = elink.call(S.fix,'/home/philbobsqp/Work/julia')$tree

# convert to ape phylogeny for plotting
N=length(mixedlabels)
Zp <- as.phylogg(Z,N+num.sources,c(rep('X',num.sources),rep('o',N)))

plot.phylo(reorder(Zp, order = "c"),edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(1:num.sources,mixedlabels,0),type='f')

# get average co-assignment from sets of individuals in the tree - add number of baseline sources as they are in the tree
Pc = get_Pc(which(mixedlabels==5)+4,which(mixedlabels==6)+4,Zp)

# Pc with a source
Pc = get_Pc(1,which(mixedlabels==1)+4,Zp)


# add circle at average co-assignment
symbols(0,0,circles=Pc,col=colors()[30],inches=F,add=T,lwd=2)

hc = hclust(dist(data))
hcp = as.phylo(hc)
hcp$tip.label=rep('o',N)

plot.phylo(reorder(hcp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1)),tip.color=c(1:num.sources,mixedlabels),type='f')

pdf('fix simulation.pdf',width=5, height=5)
par(mar=c(5,4,2,1)+0.1)
plot(baseline_t[,1],baseline_t[,2],col=cbPalette[baselabels],pch=18,xlab='LD1',ylab='LD2')
points(mixed_t[,1],mixed_t[,2],col=cbPalette[mixedlabels],pch=17)
dev.off()

pdf('fix example.pdf',width=6, height=3,colormodel='cmyk')
par(mfrow=c(1,2))
par(mar=c(5,4,2,0)+0.1)
hist(classes[keeps,1],bins,col='grey',xlab='number of sources',main='',freq=F)
par(mar=c(1,1,1,1)+0.1)
plot.phylo(reorder(Zp, order = "c"),edge.width=1.5,cex=1.5,edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(cbPalette[c(1:num.sources,mixedlabels)],0),type='f')
dev.off()

# lda classification succes
baselda <- lda(baseline,baselabels)
mean(apply(predict(baselda,mixed)$posterior,1,which.max)==mixedlabels)

# compare against mixture models

# conditional analysis
cond.output =MM.call(datas=mixed,conditional=T,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels,path.to.julia='/home/philbobsqp/Work/julia')

cond.class.id = as.data.frame(cond.output$class_id)
cc.fix=cond.class.id[,keeps]
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
uc.fix=ucond.class.id[,keeps]
for (i in num.sources:1){
uc.fix= rbind(rep(i,niter-burnin),uc.fix)}

Z = elink.call(uc.fix,path.to.julia='/home/philbobsqp/Work/julia')$tree
N=length(mixedlabels)
Zp <- as.phylogg(Z,N+num.sources,c(rep('X',num.sources),rep('o',N)))


plot.phylo(reorder(Zp, order = "c")
           ,edge.color=c(rep(1,length(Zp$edge.length)-1),0)
           ,tip.color=c(1:num.sources,mixedlabels,0)
           ,type='f')
