function [theta,predlik] = class_lik(iss,assign,priors,training_data,BN,sweep,theta,k_update)

% This code samples from the posterior of the Dirichlet Process Mixture
% model with conjugate Normal-Inv-Wishart prior, using the Split-Merge
% sampler proposed by Jain & Neal 2004 with improvements by Dahl 2005.
% Credit goes to Frank Wood (Columbia U.), on who's Gibbs sampler the Gibbs
% portion of the present code is based. The Split-Merge algo is my doing; no
% guarantee is given that code is correct or will produce correct results.
% Note, the code is quite efficient with sequential updating of sufficient
% stats at both levels (except for some places were I couldn't be asked to
% figure it out) - specifying GRAPHICS will take a LOT mote time, I therefore
% recommend to use real-time plots only initially to get things right,
% and then to update without plots. Use at your own discretion.

% Inputs are the data (training_data), the number of iterations (sweeps),
% the hyperparameters for the alpha parameter gamma prior (a_0,b_0), a
% prior guess at the mean vector (mu_0), k_0 measures the confidence in
% mu_0 (k_0 = 1 encodes little prior info); the initial degrees of freedom
% of the inv_Wishart prior - the higher, the more confident we are in our
% prior guess for the covarinace matrix Sigma (lambda_0)(needs to be number
% of dimensions of the data +1). To express a
% prior belief in Sigma, note that the expectation E of Sigma in this
% parametrization is (1/(v_0-p-1)*lambda_0^-1)^-1 , thus, to encode the prior set
% lambda_0 to (1/(v_0-p-1)*prior_E^-1)^-1. max_class_id is set to keep manageable
% data-structures. Set to a size well above realized samples - or at
% length(data).

% the code returns class_ids (a matrix of N*sweeps assignmets to classes at
% each iteration; be aware that label switching makes the numerical values of these variables
% arbitrary); K_record (a 1*sweeps vector of number of classes at each
% iteration); lp_record (a vector of log-probabilities proportional to the P
% of the model given the data); and alpha_record (the recorded alpha values
% at each iteration)

% Copyright Phil Neubauer@ Victoria University of Wellington, NZ, permission
% is granted to use, reproduce and modify the code at will. Please cite the
% corresponding paper(s) when using this code in any way.

% set alpha gamma prior parameters (vague, following Navarro)

a_0 = priors.a_0;
b_0 = priors.b_0;

% set normal inverse wishart hyper parameters

if sweep==1
    alpha = 1;
    mu_0 = priors.mu_0;
    k_0 = priors.k_0;
else
    alpha = theta(iss).alpha;
    k_0=theta(iss).G_0.k_0;
    mu_0 = theta(iss).G_0.mu_0;
end

v_0=priors.v_0;
lambda_0 = priors.var;



%% real start - if sweep =1 initialize all structures

if sweep==1
    
    [D, N] = size(training_data);
    theta(iss).D = D;
    theta(iss).N =N;
    
    max_class_id = N;
    % specify a memory efficient type for class_ids
    
    
    K_plus = 1;
    class_id = zeros(N,1);
    class_id(1) = 1;
    
    counts = zeros(max_class_id,1,'uint32');
    counts(1)=1;
    
    yyT = zeros(D,D,N);
    %for each ind allocate it to the closest mean and update the SSs
    
    for i=1:N
        y = training_data(:,i);
        yyT(:,:,i) = y*y';
    end
    theta(iss).yyT=yyT;
    %% else just get it from theta(iss)
else
    
    D = theta(iss).D;
    N = theta(iss).N;
    
    
    K_plus = theta(iss).K_record(sweep-1);
    
    counts = theta(iss).counts ;
    
    yyT = theta(iss).yyT;
    
    class_id=theta(iss).class_id;
    
end
d2 = D/2;
predlik = zeros(N-BN,1);

% run the sampler


%% update G_0


muu = zeros(D,K_plus);sums=0;
invsig = zeros(D,D,K_plus);
sig = zeros(D,D,K_plus);sumsig=0;
for k =1:K_plus
    
    m_Y = mean(training_data(:,class_id==k),2);
    
    n=double(counts(k));
    mu_n = k_0/(k_0+n)*mu_0 + n/(k_0+n)*m_Y;
    
    sum_s = sum(yyT(:,:,class_id==k),3);
    SS = (sum_s - n*(m_Y*m_Y'));
    zm_Y = m_Y-mu_0;
    lambda_n = lambda_0 + SS + k_0*n/(k_0+n)*(zm_Y)*(zm_Y)';
    
    v_n=v_0+n;
    sig(:,:,k)=iwishrnd(lambda_n,v_n);
    
    muu(:,k)=mvnrnd(mu_n,sig(:,:,k)/(k_0+n));
    
    invsig(:,:,k)=inv(sig(:,:,k));
    
    sums=sums+invsig(:,:,k)*muu(:,k);
    sumsig=sumsig+invsig(:,:,k);
    
end

meansig=inv(sumsig);

mu_0 = randnorm(1,sumsig\sums,[],meansig./k_0);

if k_update
    sums=0;
    for k=1:K_plus
        sums = sums+ (muu(:,k)-mu_0)'*invsig(:,:,k)*(muu(:,k)-mu_0);
    end
    
    k_0 = randgamma((1+double(K_plus))/2)/((1/(1+sums))/2);
    
end
    
%% p(y  under prior andalone)=(\int F(y|\theta(iss))dG_0(\theta(iss))

p_under_prior_alone = zeros(N,1);

pc_max_ind = N;
pc_gammaln_by_2 = 1:pc_max_ind;
pc_gammaln_by_2 = gammaln(pc_gammaln_by_2/2);
pc_log = reallog(1:pc_max_ind);


Sigma = (lambda_0*(k_0+1)/(k_0*(v_0-D+1)))';
v = v_0-D+1;
mu = mu_0;
log_det_Sigma = log(det(Sigma));
inv_Sigma = Sigma^-1;
vd = v+D;

for i=1:N
    y = training_data(:,i);
    
    lp = pc_gammaln_by_2(vd) - (pc_gammaln_by_2(v) + d2*pc_log(v) + ...
        d2*reallog(pi)) - .5*log_det_Sigma-...
        (vd/2)*reallog(1+(1/v)*(y-mu)'*inv_Sigma*(y-mu));
    
    p_under_prior_alone(i) = lp;
    
end

%% Gibbs sample Radford Neal's algo 2 for conjugate priors
cnt=0;
for i=randsample(N,N)'
    
    if i==1
        continue
    end
    
    y = training_data(:,i);
    
    old_class_id= class_id(i);
    
    if old_class_id ~= 0
        counts(old_class_id) = counts(old_class_id) -1;
        
        if counts(old_class_id)==0
            % delete the source and compact all data structures
            
            hits = class_id>=old_class_id;
            class_id(hits) = class_id(hits)-1;
            K_plus = K_plus-1;
            
            hits = [1:old_class_id-1 old_class_id+1:(K_plus+1)];
            
            counts(1:K_plus) = counts(hits);
            counts(K_plus+1) = 0;
            
        end
    end
    
    % complete the CRP prior with new table prob.
    priornum=[];
    priornum(counts(1:K_plus)~=0) = double(counts(counts(1:K_plus)~=0));
    priornum(counts(1:K_plus)==0) = repmat(alpha,sum(counts(1:K_plus)==0),1)/(sum(counts(1:K_plus)==0)+1);
        
    prior = [priornum';alpha/(sum(counts(1:K_plus)==0)+1)]/(sum(counts)-1+alpha);
    
    likelihood = zeros(length(prior),1);
    
    % as per Radford's Alg. 3 we will compute the posterior predictive
    % probabilities in two scenerios, 1) we will evaluate the
    % likelihood of sitting at all of the existing tables by computing
    % the probability of the datapoint under the posterior predictive
    % distribution with all points sitting at that table considered and
    % 2) we will compute the likelihood of the point under the
    % posterior predictive distribution with no observations
    
    for ell = 1:K_plus
        
        % the log likelihood for class ell
        likelihood(ell) = log(mvnpdf(y,muu(:,ell),sig(:,:,ell)));
        
    end
    likelihood(K_plus+1) = p_under_prior_alone(i);
    
    likelihoods = exp(likelihood-max(likelihood));
    likelihoods = likelihoods/sum(likelihoods);
    
    % compute the posterior over seating assignment for datum i
    posterior = prior.*likelihoods; % this is actually a proportionality
    % normalize the posterior
    posterior = posterior/sum(posterior);
    
    % pick the new table
    cdf = cumsum(posterior);
    rn = rand;
    
    new_class_id = find((cdf>rn)==1,1);
    
    if assign(i)==1
        class_id(i) = new_class_id;
        counts(new_class_id) = counts(new_class_id)+1;
    end
    
    %newcounts(new_class_id) = newcounts(new_class_id)+1;
    
    if class_id(i) == K_plus+1
        K_plus = K_plus+1;
        
        mu_n = k_0/(k_0+1)*mu_0 + n/(k_0+1)*y;
        
        zm_Y = y-mu_0;
        lambda_n = lambda_0 +  k_0*n/(k_0+1)*(zm_Y)*(zm_Y)';
        
        v_n=v_0+1;
        sig(:,:,K_plus)=iwishrnd(lambda_n,v_n);
        
        muu(:,K_plus)=mvnrnd(mu_n,sig(:,:,k)/(k_0+1));
        
    end
    
    if i>BN
        cnt=cnt+1;
        predlik(cnt)=likelihood(new_class_id);
    end
    
end

%% MCMC Alpha update code

  nu = betarnd(alpha+1,N);
    pis=(a_0+K_plus-1)/(a_0+K_plus-1+N*(b_0-log(nu)));
    alpha = (pis*(randgamma(a_0+K_plus)/(b_0-log(nu))))+((1-pis)*(randgamma(a_0+K_plus-1)/(b_0-log(nu))));
    
   
%% record the current parameters values

theta(iss).K_record(sweep) = K_plus;
theta(iss).alpha = alpha;
theta(iss).counts = counts;
theta(iss).G_0.k_0 = k_0;
theta(iss).G_0.mu_0 = mu_0;
theta(iss).class_id=class_id;



end