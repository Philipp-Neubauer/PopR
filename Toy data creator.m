%%%%%%%%%%% Toy data example %%%%%%%%%%%%

% DPM toolbox by P.Neubauer -  distributed under GNU GPL licence V3.

% this toolbox works with the lightspeed toolbox by Tom Minka,
% which can be found here:
% http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
% please install according to instructions before attempting to run any
% code

% provide your path to lightspeed here
path_to_LS = 'C:\Users\Philipp\Work\Otolith_Code\Bayesian_NP\lightspeed\lightspeed';

% load lightspeed
addpath(path_to_LS)

%% task 1 - estiamte correct number of sources without a baseline


%  multi-normal source distributions - in an ideal world...

num_sources = 3; % 'true' number of sources
num_elements = 5; % number of elements
num_per_source = 50; % individuals per source

sep = 4; % separation of means
means = randnorm(num_sources,repmat(0,num_elements,1),diag(repmat(sep,num_elements,1)));

data=zeros(num_per_source*num_sources,num_elements);
label=[];
a=0;
for i=1:num_sources
    data((a*num_per_source+1):(a*num_per_source+num_per_source),:)=randnorm(num_per_source,means(:,i))';
    label((a*num_per_source+1):(a*num_per_source+num_per_source))=i;
    a=a+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now use the DPM %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

data_DPM = data-repmat(mean(data),num_per_source*num_sources,1);

% priors is a structure with 6 elements

% prior for gamma
priors.a_0  = 1;
priors.b_0  = 1;

%Inv-Wishart prior

% degrees of freedom - at least num_elements+1
% higher v_0 equates to a more informative prior
priors.v_0  = num_elements+1;

% prior covariance - try a number of options - this shows how sensitive the
% model/data combination is to prior choice...

% play with the variance to see how the variance and source separation
% interact (use same var for each element.) - for good source separation
% this should not be very important i.e. the prior should not influence the
% number of sources. This changes when sources are not easily identifyable

vars = 1; % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
covva = diag(repmat(vars,num_elements,1));

priors.var  = inv((1/(priors.v_0-num_elements))*inv(1*covva));
% certainty about the mean...keep it low in the example
priors.k_0  = 0.01;
% prior mean
priors.mu_0 = mean(data_DPM)';

% data must be in p*n (not n*p as simulated)

niter=1000; % number of MCMC iterations - choose so that the second significant digit of E(K^+)stabilizes
intro=5; % number of Gibbs scans before starting the split merge - no need to tweak
plots=0; % plotting on ? set to 1 for yes, 0 for no


[class_id, classes,knots, gammas] = split_merge_sampler(data',priors,niter,intro,1,1,plots,path_to_LS);



% display the histogram of the number of sources

burnin = 100; % number of iterations to discard


hist(classes(burnin:niter),min(classes):bin_width:max(classes),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

% now create the exact linkage tree and display

S=[class_id(:,burnin:niter)];
Z = elink(S);
plotmod(phytree(Z))


%%%% playing with the covariance shows the sensitivity to the prior
%%%% component variance; one needs to define a source in terms of variance
%%%% about a mean to be able to estimate its number. I would personally not
%%%% recommend trying to estimate the number of sources without any prior information as it must fail (...the
%%%% sum of two normals is a normal...or inversely, the normal distribution is infinitiely divisible...)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Using the DPM with a baseline %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Case 1 - a fixed, multi-NORMAL baseline %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% same as before, set up the data
num_sources = 5 ; % sampled number of sources,
extra_sources = 3; % extra baseline
num_elements = 4; % number of elements
num_per_source(1:num_sources) = repmat(50,num_sources,1); % number of samples per source
num_per_source(num_sources+1:num_sources+extra_sources) = randi(30,extra_sources,1)+5; % number of samples per extra_sources

sep = 3; % separation of means
means = randnorm(num_sources+extra_sources,repmat(0,num_elements,1),diag(repmat(sep,num_elements,1)));

data=zeros(sum(num_per_source),num_elements);
label=[];
a=0;
for i=1:(num_sources+extra_sources)
    a=a+1;
    data(a:(a+num_per_source(i)-1),:)=randnorm(num_per_source(i),means(:,i))';
    label(a:(a+num_per_source(i))-1)=i;
    a=a+num_per_source(i)-1;
end

% visualize naively using PCA

[~,scores]=princomp(data);

% how many sources are there?
scatter(scores(:,1),scores(:,2))

% not as easy as it is with colors...
scatter(scores(:,1),scores(:,2),[],label,'LineWidth',2)

% crosses denote extra-baseline samples
hold on
for i=1:extra_sources
    scatter(scores(label==num_sources+i,1),scores(label==num_sources+i,2),'kx')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now use the DPM %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline = data(1:sum(num_per_source(1:num_sources)),:);

% add 1/4th of the samples from the baseline to the mixed set

ixs = randi(sum(num_per_source(1:num_sources)),round(sum(num_per_source(1:num_sources))/4),1);

mixed = [baseline(ixs,:); data(sum(num_per_source(1:num_sources))+1:end,:)];
baseline(ixs,:)=[];

bix = 1:sum(num_per_source(1:num_sources));
bix(ixs)=[];

mixedlabels=label([ixs' sum(num_per_source(1:num_sources))+1:end]);

mixed =  mixed-repmat(mean(baseline),size(mixed,1),1);
baseline = baseline-repmat(mean(baseline),size(baseline,1),1);

% priors is a structure with 6 elements

% prior for gamma
priors.a_0  = 1;
priors.b_0  = 1;

%Inv-Wishart prior

% degrees of freedom - at least num_elements+1
% higher v_0 equates to a more informative prior
priors.v_0  = num_elements+1;

% prior covariance - try a number of options - this shows how sensitive the
% model/data combination is to prior choice...

vari = 1; % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
priors.var  = diag(repmat(vari,num_elements,1));

% certainty about the mean...keep it low in the example
priors.k_0  = 1;
% prior mean
priors.mu_0 = mean(baseline)';

% data must be in p*n (not n*p as simulated)

niter=2000; % number of MCMC iterations - choose so that the second significant digit of E(K^+)stabilizes - usually takes at least 1000 iterations
intro=10; % number of Gibbs scans before starting the split merge and hyperparameter update - no need to tweak
plots=1; % plotting on ? set to 1 for yes, 0 for no

[class_id,~,k0lp,~,K_mix] = split_merge_sampler_fix(mixed',priors,baseline',label(bix),niter,intro,plots,1);



mean(K_mix(200:1000))

% display the histogram of the number of sources

burnin = 100; % number of iterations to discard

bin_width=1; % width of histogram bins

hist(K_mix(1:niter/thin-burnin),min(K_mix):bin_width:max(K_mix),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

% now create the exact linkage tree and display - for this one the first
% num_sources leafs are the baseline

S_fix=[repmat([1:num_sources]',1,niter/thin-burnin+1); class_id{3}(:,burnin:niter/thin)];
Z_fix = elink(S_fix);

lab=cell(1,length(mixedlabels)+num_sources);
for i=1:length(mixedlabels)+num_sources
    
    lab{i}=num2str(i);
    
end

Zfix_tree=phytree(Z_fix,lab);


% recolor tree to show the baseline (colored branches & leafs (bold, larger font), with nodes collapsed into
% one since the co-assignment proability of the baseline is 1) and mixed
% sample (colored leafs).

h=plotmod(Zfix_tree,'Type','radial');
LLs=cell(1,num_sources);

load('mycmap.mat');
cspace=linspace(1,size(cmap,1),(num_sources+extra_sources));
    
for i=1:num_sources+extra_sources
    
      
    if i<=num_sources
        LLs{i}=[num2str(i)];
        
        j=getbyname(Zfix_tree,LLs{i},'Exact','true');
        set(h.BranchLines(sum(j,2)>=1),'Color',cmap(cspace(i),:),'LineWidth',2)
        set(h.leafNodeLabels(sum(j,2)>=1),'Color',cmap(cspace(i),:),'FontWeight','bold','FontSize',14)
    end
    
    LL_group=cell(1,sum((mixedlabels==i)));
    
    for j=find(mixedlabels==i)
        LL_group{j}=[num2str(num_sources+j)];
    end
    
    jay=getbyname(Zfix_tree,LL_group,'Exact','true');
    
    set(h.leafNodeLabels(sum(jay,2)>=1),'Color',cmap(cspace(i),:))
    
end

% compare against PCA - should be the same colors, diamonds are the mixed
% sample, circles are the baseline,

bl_scores=scores(bix,:);
ms_scores=scores([ixs' sum(num_per_source(1:num_sources))+1:end],:);

figure

scatter(bl_scores(:,1),bl_scores(:,2),[],label(bix),'LineWidth',2)
hold on
scatter(ms_scores(:,1),ms_scores(:,2),[],mixedlabels','d','LineWidth',2)

% crosses denote extra-baseline samples

for i=1:extra_sources
    scatter(ms_scores(mixedlabels==num_sources+i,1),ms_scores(mixedlabels==num_sources+i,2),'kx','LineWidth',2)
end

set(gcf,'Colormap',cmap)

ylabel('PC 2')
xlabel('PC 1')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Case 2 - a fixed, but multi-MODAL baseline %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% same as before, set up the data
num_sources = 5 ; % sampled number of sources,
extra_sources = 3; % extra baseline
num_elements = 4; % number of elements
num_per_source(1:num_sources) = repmat(50,num_sources,1); % number of samples per source
num_per_source(num_sources+1:num_sources+extra_sources) = randi(30,extra_sources,1)+5; % number of samples per extra_sources

sep = 3; % separation of source-means
i_sep = 2; % intra-source separation of means
means = randnorm(num_sources+extra_sources,repmat(0,num_elements,1),diag(repmat(sep,num_elements,1)));

data=zeros(sum(num_per_source),num_elements);
label=[];fabel=[];
a=0;b=0;c=0;
for i=1:(num_sources+extra_sources)
    a=a+1;
    n_comp(i)=randi(2+1);
    ixxs = randi(num_per_source(i),n_comp(i)-1,1);
    ixxs=[sort(ixxs)' num_per_source(i)];
    % draw a random number of means for that source
    for j=1:n_comp(i)
        ixxss(j)=ixxs(j);
        if j>1
            ixxss(j)=ixxs(j)-ixxs(j-1);
        end
        b=b+1;c=c+1;
        data(b:(b+ixxss(j)-1),:)=randnorm(length(b:(b+ixxss(j)-1)),randnorm(1,means(:,i),i_sep),i_sep/2)';
        fabel(b:(b+ixxss(j)-1)) = c;
        b=b+ixxss(j)-1;
    end
    label(a:(a+num_per_source(i))-1)=i;
    
    a=a+num_per_source(i)-1;
end

% visualize naively using PCA

[~,scores]=princomp(data);

% how many sources are there ?
%scatter(scores(:,1),scores(:,2))

% not as easy as it is with colors...% crosses denote extra-baseline samples
subplot(2,1,1)
scatter(scores(:,1),scores(:,2),[],fabel,'LineWidth',2)
hold on
for i=1:extra_sources
    scatter(scores(label==num_sources+i,1),scores(label==num_sources+i,2),'kx')
end
subplot(2,1,2)
scatter(scores(:,1),scores(:,2),[],label,'LineWidth',2)
hold on
for i=1:extra_sources
    scatter(scores(label==num_sources+i,1),scores(label==num_sources+i,2),'kx')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now use the DPM %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline = data(1:sum(num_per_source(1:num_sources)),:);

% add 1/4th of the samples from the baseline to the mixed set

ixs = randi(sum(num_per_source(1:num_sources)),round(sum(num_per_source(1:num_sources))/4),1);

mixed = [baseline(ixs,:); data(sum(num_per_source(1:num_sources))+1:end,:)];
baseline(ixs,:)=[];

bix = 1:sum(num_per_source(1:num_sources));
bix(ixs)=[];

mixedlabels=label([ixs' sum(num_per_source(1:num_sources))+1:end]);

mixed =  mixed-repmat(mean(baseline),size(mixed,1),1);
baseline = baseline-repmat(mean(baseline),size(baseline,1),1);

% priors is a structure with 6 elements

% prior for gamma
priors.a_0  = 1;
priors.b_0  = 1;

%Inv-Wishart prior

% degrees of freedom - at least num_elements+1
% higher v_0 equates to a more informative prior
priors.v_0  = num_elements+1;

% prior covariance - try a number of options - this shows how sensitive the
% model/data combination is to prior choice...

var = 1; % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
covva = diag(repmat(var,num_elements,1));

priors.var  = inv((1/(priors.v_0-num_elements))*inv(1*covva));

% certainty about the mean...this is just an initial value since it is
% estimated
priors.k_0  = 1;
% prior mean
priors.mu_0 = mean(baseline)';

% data must be in p*n (not n*p as simulated)

niter=3000; % number of MCMC iterations - choose so that the second significant digit of E(K^+)stabilizes
intro=20; % number of Gibbs scans before starting the split merge - no need to tweak
plots=1; % plotting on ? set to 1 for yes, 0 for no

[class_id, K_record, lp_record, alpha_record,K_mix] = split_merge_sampler_fix(mixed',priors,baseline',label(bix),niter,intro,plots);

mean(K_mix(200:1000))

% display the histogram of the number of sources

burnin = 100; % number of iterations to discard

bin_width=1; % width of hist bins

hist(K_mix(burnin:niter),min(K_mix):bin_width:max(K_mix),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

% now create the exact linkage tree and display - for this one the first
% num_sources leafs are the baseline

S_fix=[repmat([1:num_sources]',1,niter-burnin+1); class_id(:,burnin:niter)];
Z_fix = elink(S_fix);

% create leaf labels
lab=cell(1,length(mixedlabels)+num_sources);
for i=1:length(mixedlabels)+num_sources
    
    lab{i}=num2str(i);
    
end

Z_tree=phytree(Z_fix,lab);

% recolor tree to show the baseline (colored branches & leafs, with nodes collapsed into
% one since the co-assignment proability of the baseline is 1) and mixed
% sample (colored leafs) that are from 'sampled' baseline locations.

h=plotmod(Z_tree,'Type','radial');
LLs=cell(1,num_sources);
load('mycmap.mat')
cspace=linspace(1,size(cmap,1),(num_sources+extra_sources));

for i=1:(num_sources+extra_sources)
    
       
    
    if i<=num_sources
        LLs{i}=[num2str(i)];
        
        j=getbyname(Z_tree,LLs{i},'Exact','true');
        set(h.BranchLines(sum(j,2)>=1),'Color',cmap(cspace(i),:),'LineWidth',2)
        set(h.leafNodeLabels(sum(j,2)>=1),'Color',cmap(cspace(i),:),'FontWeight','bold','FontSize',14)
    end
    
    LL_group=cell(1,sum((mixedlabels==i)));
    
    for j=find(mixedlabels==i)
        LL_group{j}=[num2str(num_sources+j)];
    end
    
    jay=getbyname(Z_tree,LL_group,'Exact','true');
    
    set(h.leafNodeLabels(sum(jay,2)>=1),'Color',cmap(cspace(i),:))
    
end

% compare against PCA - should be the same colors, diamonds are the mixed
% sample, circles are the baseline,

bl_scores=scores(bix,:);
ms_scores=scores([ixs' sum(num_per_source(1:num_sources))+1:end],:);

figure
scatter(bl_scores(:,1),bl_scores(:,2),[],label(bix),'LineWidth',2)
hold on
scatter(ms_scores(:,1),ms_scores(:,2),[],mixedlabels','d','LineWidth',2)

% crosses denote extra-baseline samples

for i=1:extra_sources
    scatter(ms_scores(mixedlabels==num_sources+i,1),ms_scores(mixedlabels==num_sources+i,2),'kx','LineWidth',2)
end
ylabel('PC 2')
xlabel('PC 1')
set(gcf,'Colormap',cmap)


% different representation, shows multimodality in the tree - same as Huelsenbck and Andolfatto 2006
Distances= 1-pdist(Z_tree,'nodes','leaves','squareform',true)/2;
imagesc(Distances)
axis off

colorbar
hold on
% plotting index
j=getbyname(Z_tree,{lab{1:num_sources}},'Exact','true');
ix = find(sum(j,2)>=1);
for i=1:length(find(sum(j,2)>=1))
    
    plot([ix(i)-0.5 ix(i)-0.5],[ix(i)-0.5 ix(i)+0.5],'color','w','LineWidth',5)
    plot([ix(i)+0.5 ix(i)+0.5],[ix(i)-0.5 ix(i)+0.5],'color','w','LineWidth',5)
    
    plot([ix(i)-0.5 ix(i)+0.5],[ix(i)-0.5 ix(i)-0.5],'color','w','LineWidth',5)
    plot([ix(i)-0.5 ix(i)+0.5],[ix(i)+0.5 ix(i)+0.5],'color','w','LineWidth',5)
    
end



