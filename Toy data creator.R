########### Toy data example ############

# DPM toolbox by P.Neubauer -  distributed under GNU GPL licence V3.

require(MASS)
require(ape)
source("Julia_call_function.R")
source("elink.call.R")
source("convert_Z_to_phylo.R")

## task 1 - estiamte correct number of sources without a baseline

#  multi-normal source distributions - in an ideal world...

num.sources = 3  # 'true' number of sources
num.elements = 5 # number of elements
num.per.source = 30 # individuals per source

sep =10 # separation of means
means = mvrnorm(num.sources,rep(0,num.elements),diag(rep(sep,num.elements)))

data=matrix(NA,num.per.source*num.sources,num.elements)
label=rep(NA,num.per.source*num.sources)
a=0
for (i in 1:num.sources){
    data[(a*num.per.source+1):(a*num.per.source+num.per.source),] <-mvrnorm(num.per.source,means[i,],diag(rep(1,num.elements)))
    label[(a*num.per.source+1):(a*num.per.source+num.per.source)] <- i
    a=a+1
}

scores = princomp(data)$scores

pdf('./Plots/very easy example.pdf',colormodel='cmyk')
plot(scores[,1],scores[,2],col=label,xlab='PCA1',ylab='PCA2')
dev.off()
##########################
### now use the DPM ######
##########################

data.DPM = data-colMeans(data)

# prior for gamma - a gamma(1,1) is reasonably broad,but wider priors usually don't make much of a difference
a.0  = 0.1
b.0  = 0.1

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

vars = 1 # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
# adjust prior by degrees of freedom to get the right expected value
lambda.0 = diag(rep(vars*(v.0-num.elements),num.elements))

# this bit does the 'adaptive' learning of 'reasonable' prior covariance.
# first run the model using, naively, the cov of all the data
lambda.0=var(data.DPM)
# skip this part for the first run, then chack which source has the lowest determinant (this is an arbitrary criterion) and set lambda.0 to the covariance of that class
for (k in sort(unique(apply(class.id[,burnin:niter],1,median)))){
cat(det(cov(data.DPM[apply(class.id[,burnin:niter],1,median)==k,])),'\n')}

# set the number after the == to the source with the lowest determinant
lambda.0=var(data.DPM[apply(class.id[,burnin:niter],1,median)==1,])



# certainty about the mean...keep it low in the example
k.0  = 1
# prior mean
mu.0 = colMeans(data.DPM)

# number of iterations per processor
num.iters=1000
# numebr of parallel processing jobs
np=1
# thinning of the marcov chain
thin=2
# total number of kept iterations
niter=np*num.iters/thin

output = DPM.call(datas=data.DPM,iters=num.iters,thin=thin,np=np, path.to.julia='/home/philbert/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)

# display the histogram of the number of sources

burnin = 100  # number of (kept!) iterations to discard
bins = (min(classes)-0.5):(max(classes)+0.5) # histogram bins

pdf('./Plots/hist very  easy example.pdf',colormodel='cmyk')
hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)
dev.off()


# now create the exact linkage tree and display

S=class.id[,burnin:niter]
Z = elink.call(S,path.to.julia='/home/philbert/julia')$tree
N=num.per.source*num.sources
Zp <- as.phylogg(Z,num.per.source*num.sources,rep('o',N))

#pdf('./Plots/tree very easy example.pdf')
plot.phylo(reorder(Zp, order = "c"),edge.color=c(rep(1,length(Zp$edge.length)-1),0),tip.color=c(label,0),type='f')
#text(0.03,0,0)
#text(-0.05,0,0.05)
#axisPhylo()
#dev.off()

pdf('./Plots/clustering very easy example.pdf')
plot(hclust(dist(data.DPM)),labels=F)
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

mixed =  mixed-colMeans(baseline)
baseline = baseline-colMeans(baseline)


# prior for gamma - a gamma(1,1) is reasonably broad,but wider priors usually don't make much of a difference
a.0  = 1
b.0  = 1

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

vars = 10 # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
lambda.0 = diag(rep(vars,num.elements))

# certainty about the mean...keep it low in the example
k.0  = 0.01
# prior mean
mu.0 = colMeans(baseline)

# number of iterations per processor
num.iters=10000
# numebr of parallel processing jobs
np=1
# thinning of the marcov chain
thin=10
# total number of kept iterations
niter=np*num.iters/thin

output =DPM.call(datas=mixed,learn=T,iters=num.iters,thin=thin,np=np,baseline=baseline,labels=baselabels)

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)

# display the histogram of the number of sources

burnin = 100  # number of (kept!) iterations to discard
bins = (min(classes)-0.5):(max(classes)+0.5) # histogram bins

pdf('./Plots/hist fix very  easy example.pdf',colormodel='cmyk')
hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)
dev.off()

# now create the exact linkage tree and display - for this one the first
# num.sources leafs are the baseline
S.fix=class.id[,burnin:niter]
for (i in num.sources:1){
S.fix= rbind(rep(i,niter-burnin),S.fix)}

Z = elink.call(S.fix,path.to.julia='/home/philbert/julia')$tree
N=length(mixedlabels)
Zp <- as.phylogg(Z,N+num.sources,c(rep('X',num.sources),rep('o',N)))

pdf('./Plots/tree fix very easy example_var10.pdf')
plot.phylo(reorder(Zp, order = "c")
           ,edge.color=c(rep(1,length(Zp$edge.length)-1),0)
           ,tip.color=c(1:num.sources,mixedlabels,0)
           ,type='f')
#text(0.03,0,0)
#text(-0.05,0,0.05)
#axisPhylo()
dev.off()

pdf('./Plots/clustering very easy example.pdf')
plot(hclust(dist(mixed)),labels=F)
dev.off()




