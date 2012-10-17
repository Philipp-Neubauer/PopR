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

sep = 6 # separation of means
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

vars = 1 # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
lambda.0 = diag(rep(vars,num.elements))

# certainty about the mean...keep it low in the example
k.0  = 0.01
# prior mean
mu.0 = colMeans(data.DPM)

# number of iterations per processor
num.iters=5000
# numebr of parallel processing jobs
np=2
# thinning of the marcov chain
thin=10
# total number of kept iterations
niter=np*num.iters/thin

output = DPM.call(baseline=F,iters=num.iters,thin=thin,np=np, path.to.julia='/home/philbert/julia')

# these are the source allocations for all kept MCMC iterations
class.id = as.data.frame(output$class_id)
# and just the number of sources per iteration
classes = as.data.frame(output$K_record)

# display the histogram of the number of sources

burnin = 10  # number of (kept!) iterations to discard
bins = (min(classes)-0.5):(max(classes)+0.5) # histogram bins

pdf('./Plots/hist very  easy example.pdf',colormodel='cmyk')
hist(classes[burnin:niter,1],bins,col='grey',xlab='number of sources',main='',freq=F)
dev.off()

# now create the exact linkage tree and display

S=class.id[,burnin:niter]
Z = elink.call(S,path.to.julia='/home/philbert/julia')$tree
N=num.per.source*num.sources
Zp <- as.phylogg(Z,num.per.source*num.sources,rep('o',N))

pdf('./Plots/tree very easy example.pdf')
plot.phylo(reorder(Zp, order = "c")
           ,edge.color=c(rep(1,length(Zp$edge.length)-1),0)
           ,tip.color=c(label,0)
           ,type='f')
#text(0.03,0,0)
#text(-0.05,0,0.05)
#axisPhylo()
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
#################################################
#### Case 1 - a fixed, multi-NORMAL baseline ####
#################################################
#########################################################

# same as before, set up the data
num.sources = 4  # sampled number of sources,
extra.sources = 2 # extra baseline
num.elements = 5 # number of elements
num.per.source = rep(30,num.sources) # number of samples per source
num.per.source = c(num.per.source,sample.int(30,extra.sources)+5) # number of samples per extra.sources
as=num.sources+extra.sources

sep = 7 # separation of means
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

pdf('./Plots/easy example fix.pdf',colormodel='cmyk')
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

vars = 1 # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
lambda.0 = diag(rep(vars,num.elements))

# certainty about the mean...keep it low in the example
k.0  = 0.01
# prior mean
mu.0 = colMeans(data.DPM)

# number of iterations per processor
num.iters=500
# numebr of parallel processing jobs
np=2
# thinning of the marcov chain
thin=10
# total number of kept iterations
niter=np*num.iters/thin

output = DPM.call(learn=T,iters=num.iters,thin=thin,np=np, path.to.julia='/home/philbert/julia',baseline=baseline,labels=baselabels)


mean(K.mix(200:1000))

# display the histogram of the number of sources

burnin = 100 # number of iterations to discard

bin.width=1 # width of histogram bins

hist(K.mix(1:niter/thin-burnin),min(K.mix):bin.width:max(K.mix),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

# now create the exact linkage tree and display - for this one the first
# num.sources leafs are the baseline

S.fix=[repmat([1:num.sources]',1,niter/thin-burnin+1) class.id{3}(:,burnin:niter/thin)]
Z.fix = elink(S.fix)

lab=cell(1,length(mixedlabels)+num.sources)
for i=1:length(mixedlabels)+num.sources
    
    lab{i}=num2str(i)
    
end

Zfix.tree=phytree(Z.fix,lab)


# recolor tree to show the baseline (colored branches & leafs (bold, larger font), with nodes collapsed into
# one since the co-assignment proability of the baseline is 1) and mixed
# sample (colored leafs).

h=plotmod(Zfix.tree,'Type','radial')
LLs=cell(1,num.sources)

load('mycmap.mat')
cspace=linspace(1,size(cmap,1),(num.sources+extra.sources))
    
for i=1:num.sources+extra.sources
    
      
    if i<=num.sources
        LLs{i}=[num2str(i)]
        
        j=getbyname(Zfix.tree,LLs{i},'Exact','true')
        set(h.BranchLines(sum(j,2)>=1),'Color',cmap(cspace(i),:),'LineWidth',2)
        set(h.leafNodeLabels(sum(j,2)>=1),'Color',cmap(cspace(i),:),'FontWeight','bold','FontSize',14)
    end
    
    LL.group=cell(1,sum((mixedlabels==i)))
    
    for j=find(mixedlabels==i)
        LL.group{j}=[num2str(num.sources+j)]
    end
    
    jay=getbyname(Zfix.tree,LL.group,'Exact','true')
    
    set(h.leafNodeLabels(sum(jay,2)>=1),'Color',cmap(cspace(i),:))
    
end

# compare against PCA - should be the same colors, diamonds are the mixed
# sample, circles are the baseline,

bl.scores=scores(bix,:)
ms.scores=scores([ixs' sum(num.per.source(1:num.sources))+1:end],:)

figure

scatter(bl.scores(:,1),bl.scores(:,2),[],label(bix),'LineWidth',2)
hold on
scatter(ms.scores(:,1),ms.scores(:,2),[],mixedlabels','d','LineWidth',2)

# crosses denote extra-baseline samples

for i=1:extra.sources
    scatter(ms.scores(mixedlabels==num.sources+i,1),ms.scores(mixedlabels==num.sources+i,2),'kx','LineWidth',2)
end

set(gcf,'Colormap',cmap)

ylabel('PC 2')
xlabel('PC 1')

##
######################################################
#### Case 2 - a fixed, but multi-MODAL baseline #####
#####################################################

# same as before, set up the data
num.sources = 5  # sampled number of sources,
extra.sources = 3 # extra baseline
num.elements = 4 # number of elements
num.per.source(1:num.sources) = repmat(50,num.sources,1) # number of samples per source
num.per.source(num.sources+1:num.sources+extra.sources) = randi(30,extra.sources,1)+5 # number of samples per extra.sources

sep = 3 # separation of source-means
i.sep = 2 # intra-source separation of means
means = randnorm(num.sources+extra.sources,repmat(0,num.elements,1),diag(repmat(sep,num.elements,1)))

data=zeros(sum(num.per.source),num.elements)
label=[]fabel=[]
a=0b=0c=0
for i=1:(num.sources+extra.sources)
    a=a+1
    n.comp(i)=randi(2+1)
    ixxs = randi(num.per.source(i),n.comp(i)-1,1)
    ixxs=[sort(ixxs)' num.per.source(i)]
    # draw a random number of means for that source
    for j=1:n.comp(i)
        ixxss(j)=ixxs(j)
        if j>1
            ixxss(j)=ixxs(j)-ixxs(j-1)
        end
        b=b+1c=c+1
        data(b:(b+ixxss(j)-1),:)=randnorm(length(b:(b+ixxss(j)-1)),randnorm(1,means(:,i),i.sep),i.sep/2)'
        fabel(b:(b+ixxss(j)-1)) = c
        b=b+ixxss(j)-1
    end
    label(a:(a+num.per.source(i))-1)=i
    
    a=a+num.per.source(i)-1
end

# visualize naively using PCA

[~,scores]=princomp(data)

# how many sources are there ?
#scatter(scores(:,1),scores(:,2))

# not as easy as it is with colors...# crosses denote extra-baseline samples
subplot(2,1,1)
scatter(scores(:,1),scores(:,2),[],fabel,'LineWidth',2)
hold on
for i=1:extra.sources
    scatter(scores(label==num.sources+i,1),scores(label==num.sources+i,2),'kx')
end
subplot(2,1,2)
scatter(scores(:,1),scores(:,2),[],label,'LineWidth',2)
hold on
for i=1:extra.sources
    scatter(scores(label==num.sources+i,1),scores(label==num.sources+i,2),'kx')
end



##########################
### now use the DPM ######
##########################

baseline = data(1:sum(num.per.source(1:num.sources)),:)

# add 1/4th of the samples from the baseline to the mixed set

ixs = randi(sum(num.per.source(1:num.sources)),round(sum(num.per.source(1:num.sources))/4),1)

mixed = [baseline(ixs,:) data(sum(num.per.source(1:num.sources))+1:end,:)]
baseline(ixs,:)=[]

bix = 1:sum(num.per.source(1:num.sources))
bix(ixs)=[]

mixedlabels=label([ixs' sum(num.per.source(1:num.sources))+1:end])

mixed =  mixed-repmat(mean(baseline),size(mixed,1),1)
baseline = baseline-repmat(mean(baseline),size(baseline,1),1)

# priors is a structure with 6 elements

# prior for gamma
priors.a.0  = 1
priors.b.0  = 1

#Inv-Wishart prior

# degrees of freedom - at least num.elements+1
# higher v.0 equates to a more informative prior
priors.v.0  = num.elements+1

# prior covariance - try a number of options - this shows how sensitive the
# model/data combination is to prior choice...

var = 1 # consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
covva = diag(repmat(var,num.elements,1))

priors.var  = inv((1/(priors.v.0-num.elements))*inv(1*covva))

# certainty about the mean...this is just an initial value since it is
# estimated
priors.k.0  = 1
# prior mean
priors.mu.0 = mean(baseline)'

# data must be in p*n (not n*p as simulated)

niter=3000 # number of MCMC iterations - choose so that the second significant digit of E(K^+)stabilizes
intro=20 # number of Gibbs scans before starting the split merge - no need to tweak
plots=1 # plotting on ? set to 1 for yes, 0 for no

[class.id, K.record, lp.record, alpha.record,K.mix] = split.merge.sampler.fix(mixed',priors,baseline',label(bix),niter,intro,plots)

mean(K.mix(200:1000))

# display the histogram of the number of sources

burnin = 100 # number of iterations to discard

bin.width=1 # width of hist bins

hist(K.mix(burnin:niter),min(K.mix):bin.width:max(K.mix),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

# now create the exact linkage tree and display - for this one the first
# num.sources leafs are the baseline

S.fix=[repmat([1:num.sources]',1,niter-burnin+1) class.id(:,burnin:niter)]
Z.fix = elink(S.fix)

# create leaf labels
lab=cell(1,length(mixedlabels)+num.sources)
for i=1:length(mixedlabels)+num.sources
    
    lab{i}=num2str(i)
    
end

Z.tree=phytree(Z.fix,lab)

# recolor tree to show the baseline (colored branches & leafs, with nodes collapsed into
# one since the co-assignment proability of the baseline is 1) and mixed
# sample (colored leafs) that are from 'sampled' baseline locations.

h=plotmod(Z.tree,'Type','radial')
LLs=cell(1,num.sources)
load('mycmap.mat')
cspace=linspace(1,size(cmap,1),(num.sources+extra.sources))

for i=1:(num.sources+extra.sources)
    
       
    
    if i<=num.sources
        LLs{i}=[num2str(i)]
        
        j=getbyname(Z.tree,LLs{i},'Exact','true')
        set(h.BranchLines(sum(j,2)>=1),'Color',cmap(cspace(i),:),'LineWidth',2)
        set(h.leafNodeLabels(sum(j,2)>=1),'Color',cmap(cspace(i),:),'FontWeight','bold','FontSize',14)
    end
    
    LL.group=cell(1,sum((mixedlabels==i)))
    
    for j=find(mixedlabels==i)
        LL.group{j}=[num2str(num.sources+j)]
    end
    
    jay=getbyname(Z.tree,LL.group,'Exact','true')
    
    set(h.leafNodeLabels(sum(jay,2)>=1),'Color',cmap(cspace(i),:))
    
end

# compare against PCA - should be the same colors, diamonds are the mixed
# sample, circles are the baseline,

bl.scores=scores(bix,:)
ms.scores=scores([ixs' sum(num.per.source(1:num.sources))+1:end],:)

figure
scatter(bl.scores(:,1),bl.scores(:,2),[],label(bix),'LineWidth',2)
hold on
scatter(ms.scores(:,1),ms.scores(:,2),[],mixedlabels','d','LineWidth',2)

# crosses denote extra-baseline samples

for i=1:extra.sources
    scatter(ms.scores(mixedlabels==num.sources+i,1),ms.scores(mixedlabels==num.sources+i,2),'kx','LineWidth',2)
end
ylabel('PC 2')
xlabel('PC 1')
set(gcf,'Colormap',cmap)


# different representation, shows multimodality in the tree - same as Huelsenbck and Andolfatto 2006
Distances= 1-pdist(Z.tree,'nodes','leaves','squareform',true)/2
imagesc(Distances)
axis off

colorbar
hold on
# plotting index
j=getbyname(Z.tree,{lab{1:num.sources}},'Exact','true')
ix = find(sum(j,2)>=1)
for i=1:length(find(sum(j,2)>=1))
    
    plot([ix(i)-0.5 ix(i)-0.5],[ix(i)-0.5 ix(i)+0.5],'color','w','LineWidth',5)
    plot([ix(i)+0.5 ix(i)+0.5],[ix(i)-0.5 ix(i)+0.5],'color','w','LineWidth',5)
    
    plot([ix(i)-0.5 ix(i)+0.5],[ix(i)-0.5 ix(i)-0.5],'color','w','LineWidth',5)
    plot([ix(i)-0.5 ix(i)+0.5],[ix(i)+0.5 ix(i)+0.5],'color','w','LineWidth',5)
    
end



