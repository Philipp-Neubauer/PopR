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

varss=[0.1 1 10];
as=0;nss=0;classee=[];knotss=[];
%  multi-normal source distributions - in an ideal world...
for k=2:0.5:4
    as=as+1;
    for t=1:5
        nss=0;
        for ns = [3 6 12]
            nss=nss+1;
            num_sources = ns; % 'true' number of sources
            num_elements = 5; % number of elements
            num_per_source = 50;
            
            sep = k; % separation of means
            means = randnorm(num_sources,repmat(0,num_elements,1),diag(repmat(sep,num_elements,1)));
            
            data=zeros(num_per_source*num_sources,num_elements);
            label=[];
            a=0;
            for i=1:num_sources
                data((a*num_per_source+1):(a*num_per_source+num_per_source),:)=randnorm(num_per_source,means(:,i))';
                label((a*num_per_source+1):(a*num_per_source+num_per_source))=i;
                a=a+1;
            end
            
            for v=1:3
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
                
                % naively use the covariance matrix of all the data...obviously a bad
                % choice because it will gegerally lead to fewer sources
                covva = cov(data);
                
                % play with the variance to see how the variance and source separation
                % interact (use same var for each element.) - for good source separation
                % this should not be very important i.e. the prior should not influence the
                % number of sources. This changes when sources are not easily identifyable
                
                var = varss(v); % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
                covva = diag(repmat(var,num_elements,1));
                
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
                
                classee(as,nss,v,t)=mode(classes(200:1000));
                knotss(as,nss,v,t)=median(knots);
                
                
            end
        end
    end
end
% display the histogram of the number of sources

burnin = 100; % number of iterations to discard

save classee.mat classee

subplot(3,1,1)
mc= mode(classee(:,1,:,:),4);
image(1:5,(varss),(reshape(mc,5,3)'))
colormap('lines')
set(gca,'YTick',0:5:10)
set(gca, 'YTickLabel',varss)
ylabel('\delta_0')
xlabel('distance')
set(gca,'XTick',1:5)
colorbar('YLim',[-0.5 3.5])

subplot(3,1,2)
mc= mode(classee(:,2,:,:),4);
image(1:5,(varss),(reshape(mc,5,3)'))
colormap('lines')
set(gca,'YTick',0:5:10)
set(gca, 'YTickLabel',varss)
ylabel('\delta_0')
xlabel('distance')
set(gca,'XTick',1:5)
colorbar('YLim',[-0.5 5.5])

subplot(3,1,3)
mc= mode(classee(:,3,:,:),4);
image(1:5,(varss),(reshape(mc,5,3)'))
set(gca,'YTick',0:5:10)
colormap('lines')
set(gca, 'YTickLabel',varss)
ylabel('\delta_0')
xlabel('distance')
set(gca,'XTick',1:5)
colorbar('YLim',[1.5 8.5])
% bin_width=1; % width of hist bins

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

%% Using the DPM with a baseline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Case 1 - a fixed, multi-normal baseline %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% same as before, set up the data
num_sources = 5 ; % sampled number of sources,
extra_sources = 3; % extra baseline
num_elements = 4; % number of elements
num_per_source(1:num_sources) = repmat(50,num_sources,1); % number of samples per source
num_per_source(num_sources+1:num_sources+extra_sources) = randi(30,extra_sources,1)+5; % number of samples per extra_sources

sep = 4; % separation of means
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
priors.a_0  = 0;
priors.b_0  = 0;

%Inv-Wishart prior

% degrees of freedom - at least num_elements+1
% higher v_0 equates to a more informative prior
priors.v_0  = num_elements+1;

% prior covariance - try a number of options - this shows how sensitive the
% model/data combination is to prior choice...
for i=2
    var = vars(i); % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
    covva = diag(repmat(var,num_elements,1));
    
    priors.var  = inv((1/(priors.v_0-num_elements))*inv(1*covva));
    % certainty about the mean...keep it low in the example
    priors.k_0  = 1;
    % prior mean
    priors.mu_0 = mean(baseline)';
    
    % data must be in p*n (not n*p as simulated)
    
    niter=10000; % number of MCMC iterations - choose so that the second significant digit of E(K^+)stabilizes - usually takes at least 1000 iterations
    intro=10; % number of Gibbs scans before starting the split merge and hyperparameter update - no need to tweak
    plots=1; % plotting on ? set to 1 for yes, 0 for no
    
    [class_id{i},~,k0hp{i},~,K_mix] = split_merge_sampler_fix(mixed',priors,baseline',label(bix),niter,intro,plots,10);
    
    KMhp{i} = K_mix(100:1000);
    
end

mean(K_mix(200:1000))

for i =1:5
    kml(i)=median(KMlp{i});
   
   sml(i,:)=quantile(KMlp{i},[0.025,0.975]);
end


for i =1:5
    km(i)=median(KMhp{i});
    sm(i,:)=quantile(KMhp{i},[0.025,0.975]);
end

errorbar(1:5,km,sm(:,1),sm(:,2),'ko-')
hold on
errorbar(1:5,kml,sml(:,1),sml(:,2),'bo-')
ylabel('number of sources')
xlabel('\delta_0')
axis([0 6 0 22]);
set(gca,'XTick',1:5)
set(gca, 'XTickLabel',vars)


% display the histogram of the number of sources

burnin = 100; % number of iterations to discard

bin_width=1; % width of hist bins

hist(KMhp{1}(1:niter/thin-burnin),min(KMhp{1}):bin_width:max(KMhp{1}),'FaceColor','grey')
set(gca,'YTick',[])
ylabel('density')
xlabel('number of sources')

% now create the exact linkage tree and display - for this one the first
% num_sources leafs are the baseline

S_fix=[repmat([1:num_sources]',1,niter/thin-burnin+1); class_id{2}(:,burnin:niter/thin)];
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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Case 4 - Classification with an incomplete %%%%%%
%%%% hierarchical  baseline                     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Imagine now that the num_sources previously sampled sources are the only
% possible source 'regions', but that sampling within these regions was
% incomplete - i.e. intraregion variability is important such that the
% baseline is multi-modal, but is missing some modes due to limited
% sampling. We now classify to one of the previous regions, assuming that
% the extra-sources originate from one of these regions.

% INPUTS:

% mixture = mixed sample
% baseline = baseline samples
% baseline_sources_high = regional source labels for the baseline
% baseline_sources_low = within region source labels 0- can be the same as
% above
% can be the same as baseline_sources_high if there is no hierarchical structure
% rest as in previous cases

% OUTPUTS:
% pis  - mixing proportions
% assignments - probabilistic assigments
% reg-assignements - assignment to regional ensembles (same as above if only one group per region)
% theta - structure of parameter values at last iteration
% pzs - posterior assignment probabilities used for mixture or discriminant assignment

% Discriminant (0) or mixture estiamtion (1) (estiamtes the mixing proportions)

mix=1;

var = 1; % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
covva = diag(repmat(var,num_elements,1));

priors.var  = inv((1/(priors.v_0-num_elements))*inv(1*covva));

niter=5000;
burnin=1000;

[pis,assignments,reg_assignments,theta,pzs] = DPMcMCMC(mix,mixed',baseline',label(bix)',fabel(bix)',priors,niter,1,path_to_LS);


[n,xout]=hist(pis);

bar(xout,n,2)

cspace=linspace(1,size(cmap,1),(num_sources+extra_sources));

h = findobj(gca,'Type','patch');
for i=1:length(h)
    set(h(i),'FaceColor',cmap(cspace(i),:))
end
set(gcf,'Colormap',cmap)

set(gca,'YTick',[])
ylabel('density')
xlabel('source proportions')

% now create the exact linkage tree and display - for this one the first
% num_sources leafs are the baseline

S_class=[repmat([1:num_sources]',1,niter-burnin+1); assignments(:,burnin:niter)];
Z_class = elink(S_class);

% create leaf labels
lab=cell(1,size(assignments,1)+num_sources);
for i=1:size(assignments,1)+num_sources
    
    lab{i}=num2str(i);
    
end

Z_tree=phytree(Z_class,lab);

% recolor tree to show the baseline (colored branches & leafs, with nodes collapsed into
% one since the co-assignment proability of the baseline is 1) and mixed
% sample (colored leafs) that are from 'sampled' baseline locations.

h=plotmod(Z_tree,'Type','radial');
LLs=cell(1,num_sources);
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

% compare to LDA

mean(mixedlabels(mixedlabels<=num_sources)==mode(assignments(mixedlabels<=num_sources,:)'))

[class,ps,pps] = classify(mixed,baseline,label(bix));
mean(class(mixedlabels<=num_sources)' == mixedlabels(mixedlabels<=num_sources))

[class'; mode(assignments')]

mpzs=mean(pzs,3); ppsi=[];pzsi=[];
for i=1:length(class)
    ppsi(i) = pps(i,class(i));
    pzsi(i) = mpzs(i,class(i));
end

addpath C:\Users\Philipp\Work\Otolith_Code\discrim

F=lda(baseline,label(bix)); %do LDA on Sites

[lambda,Cv]=cvar(F);  % lambda gives the eigenvectors, Cv gives the % variance
plot(Cv)
sum(Cv(1:3))

% plot LDA  and DPMc results side by side - rounds are mixed sample, scaled
% by the psoterior classification probability

figure(1)
subplot(1,2,1)
scatter(baseline*lambda(:,1),baseline*lambda(:,2),50,label(bix),'x','linewidth',2)
hold on
scatter(mixed*lambda(:,1),mixed*lambda(:,2),50*ppsi,class,'o','linewidth',2)
subplot(1,2,2)
scatter(baseline*lambda(:,1),baseline*lambda(:,2),50,label(bix),'x','linewidth',2)
hold on
scatter(mixed*lambda(:,1),mixed*lambda(:,2),50*pzsi,mode(assignments'),'o','linewidth',2)

figure(2)
subplot(1,2,1)
imagesc(pps)
line([0.5,5.5],[64,64],'color','w','LineWidth',2)
set(gca,'XTick',1:5)
subplot(1,2,2)
imagesc(mpzs)
line([0.5,5.5],[64,64],'color','w','LineWidth',2)
set(gca,'XTick',1:5)

%% redo classification, this time all individuals come from sampled regions, jsut from unsmapled subbpops
num_per_source=[];
num_sources = 4 ; % sampled number of sources,
num_elements = 4; % number of elements
num_per_source(1:num_sources) = repmat(150,num_sources,1); % number of samples per source

sep = 3.5; % separation of source-means
i_sep = 2; % intra-source separation of means
means = randnorm(num_sources+extra_sources,repmat(0,num_elements,1),diag(repmat(sep,num_elements,1)));

data=zeros(sum(num_per_source),num_elements);
label=[];fabel=[];
a=0;b=0;c=0;
for i=1:(num_sources)
    a=a+1;
    n_comp(i)=4;
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

% some sites within regions are not sampled
omit = [3 6 12 16];
% not as easy as it is with colors...% crosses denote extra-baseline samples
subplot(2,1,1)
scatter(scores(:,1),scores(:,2),[],fabel,'LineWidth',2)
hold on
for i=1:4
    scatter(scores(fabel==omit(i),1),scores(fabel==omit(i),2),'kx')
end

subplot(2,1,2)
scatter(scores(:,1),scores(:,2),[],label,'LineWidth',2)

mixed =[];baseline = data(1:sum(num_per_source(1:num_sources)),:);fabels=fabel;olab=[];labels=label;
for i=1:4
    mixed = [mixed; data(fabel==omit(i),:)];
    olab=[olab label(fabels==omit(i))];
    labels(fabels==omit(i))=[];
    baseline(fabels==omit(i),:)=[];
    fabels(fabels==omit(i))=[];
end

% add 1/4th of the samples from the baseline to the mixed set
lbs=length(baseline);
ixs = randi(lbs,round(lbs/5),1);

mixed = [baseline(ixs,:); mixed];
baseline(ixs,:)=[];

bix = 1:(lbs);
bix(ixs)=[];

mixedlabels=[labels(ixs) olab];

mixed =  mixed-repmat(mean(baseline),size(mixed,1),1);
baseline = baseline-repmat(mean(baseline),size(baseline,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Use DMPc %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mix=1;

var = 1; % consider changing this over orders of margnitude - e.g., 0.1,1,10 and rerun the anlysis with each
covva = diag(repmat(var,num_elements,1));

priors.var  = inv((1/(priors.v_0-num_elements))*inv(1*covva));

niter=1000;
burnin=100;

[piss,ass,reg_ass,~,pzss] = DPMcMCMC(mix,mixed',baseline',labels(bix)',fabels(bix)',priors,niter,1,path_to_LS);

[n,xout]=hist(piss);

bar(xout,n,2)

cspace=linspace(1,size(cmap,1),4);

h = findobj(gca,'Type','patch');
for i=1:length(h)
    set(h(i),'FaceColor',cmap(cspace(i),:))
end
set(gcf,'Colormap',cmap)

set(gca,'YTick',[])
ylabel('density')
xlabel('source proportions')

% now create the exact linkage tree and display - for this one the first
% num_sources leafs are the baseline

S_class=[repmat([1:num_sources]',1,niter-burnin+1); ass(:,burnin:niter)];
Z_class = elink(S_class);

% create leaf labels
lab=cell(1,size(ass,1)+num_sources);
for i=1:size(ass,1)+num_sources
    
    lab{i}=num2str(i);
    
end

Z_tree=phytree(Z_class,lab);

% recolor tree to show the baseline (colored branches & leafs, with nodes collapsed into
% one since the co-assignment proability of the baseline is 1) and mixed
% sample (colored leafs) that are from 'sampled' baseline locations.

h=plotmod(Z_tree,'Type','radial');
LLs=cell(1,num_sources);
cspace=linspace(1,size(cmap,1),(num_sources));

for i=1:(num_sources)
    
   
   
    
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
scatter(bl_scores(:,1),bl_scores(:,2),[],labels(bix),'LineWidth',2)
hold on
scatter(ms_scores(:,1),ms_scores(:,2),[],mixedlabels','d','LineWidth',2)

% crosses denote extra-baseline samples

ylabel('PC 2')
xlabel('PC 1')
set(gcf,'Colormap',cmap)

% compare to LDA

mean(mixedlabels==mode(ass'))

[classs,pss,ppss] = classify(mixed,baseline,labels(bix));
mean(classs' == mixedlabels)

mpzs=mean(pzss,3); ppsi=[];pzsi=[];
for i=1:length(classs)
    ppsi(i) = ppss(i,classs(i));
    pzsi(i) = mpzs(i,classs(i));
end

addpath C:\Users\Philipp\Work\Otolith_Code\discrim

F=lda(baseline,labels(bix)); %do LDA on Sites

[lambda,Cv]=cvar(F);  % lambda gives the eigenvectors, Cv gives the % variance
plot(Cv)
sum(Cv(1:3))

% plot LDA  and DPMc results side by side - rounds are mixed sample, scaled
% by the psoterior classification probability

figure(1)
subplot(1,2,1)
scatter(baseline*lambda(:,1),baseline*lambda(:,2),50,labels(bix),'o','linewidth',2)
hold on
scatter(mixed*lambda(:,1),mixed*lambda(:,2),50*ppsi,classs,'d','linewidth',2)
scatter(mixed(classs' ~= mixedlabels,:)*lambda(:,1),mixed(classs' ~= mixedlabels,:)*lambda(:,2),50*ppsi(classs' ~= mixedlabels),'kx','linewidth',2)

subplot(1,2,2)
scatter(baseline*lambda(:,1),baseline*lambda(:,2),50,labels(bix),'o','linewidth',2)
hold on
scatter(mixed*lambda(:,1),mixed*lambda(:,2),50*pzsi,mode(ass'),'d','linewidth',2)
scatter(mixed(classs' ~= mixedlabels,:)*lambda(:,1),mixed(classs' ~= mixedlabels,:)*lambda(:,2),50*pzsi(classs' ~= mixedlabels),'kx','linewidth',2)
