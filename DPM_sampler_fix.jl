function  DPM_sampler_fix(num_iters,thin)

single_priors = dlmread("single_priors.csv",",",Float64)
single_priors=single_priors[2:end,2]

matrix_priors = dlmread("matrix_priors.csv",",",Float64)
matrix_priors=matrix_priors[2:end,2:end]

datas =dlmread("datas.csv",",",Float64)
global datas=datas[2:end,2:end]'

baseline =dlmread("baseline.csv",",",Float64)
global baseline=baseline[2:end,2:end]'

label = dlmread("labels.csv",",",Float64)
label = label[2:end,2]
label=int(label)

    load("Julia_code_support.jl")
    load("gibbs_crp_fix.jl")
    load("split-merge_fix.jl")

    uniqs=unique(label,true)
    labels=Array(Int64,length(label))
    for i=1:length(uniqs)
        labels[label.==uniqs[i]]=i
    end
    
    a_0 = single_priors[1]
    b_0= single_priors[2]

    k_0 = single_priors[3]
    v_0 = single_priors[4]

    mu_0= single_priors[5:end]

    lambda_0 =  matrix_priors
    

    alpha = 1

    # initialize structures

    (D, N) = size(datas)
    max_class_id = N

    nit=int(num_iters/thin)
    sources=max(labels)
    k_0s=Array(Float64,nit)
    K_plus = sources
    class_id = Array(Int64,N)
    class_ids = Array(Int64,(N,nit))
    K_record = Array(Int64,nit)
    K_mix_record = Array(Int64,nit)
    alpha_record = Array(Float64,nit)
    p_under_prior_alone = Array(Float64,N);
      
    # precompute student-T posterior predictive distribution constants
    const pc_max_ind = 1e5
    const pc_gammaln_by_2 = lgamma((1:pc_max_ind)/2)
    const pc_log_pi = log(pi)
    const pc_log = log(1:pc_max_ind)
    
    means = Array(Float64,(D,max_class_id))
    orgmeans = Array(Float64,(D,max_class_id))
    sum_squares = Array(Float64,(D,D,max_class_id))
    orgsum_squares = Array(Float64,(D,D,max_class_id))
    inv_cov = Array(Float64,(D,D,max_class_id))
    log_det_cov = Array(Float64,(max_class_id))
    counts = Array(Int64,max_class_id,1)
    allcounts = Array(Int64,max_class_id,1)
    ns = Array(Int64,max_class_id);
    is=unique(labels,true);
    yyT = Array(Float64,(D,D,N))
    orginv_cov=Array(Float64,(D,D,sources))
    orglog_det_cov=Array(Float64,(D,D,sources))

    for i=1:sources
    
        means[:,i] = mean(baseline[:,labels.==is[i]],2)
        orgmeans[:,i] =  means[:,i]
        allcounts[i] = size(baseline[:,labels.==is[i]],2)
        ns[i] = size(baseline[:,labels.==is[i]],2)
        orgsum_squares[:,:,i] = (baseline[:,labels.==is[i]]-repmat(mean(baseline[:,labels.==is[i]],2),1,allcounts[i]))*(baseline[:,labels.==is[i]]-repmat(mean(baseline[:,labels.==is[i]],2),1,allcounts[i]))'
        sum_squares[:,:,i] =  orgsum_squares[:,:,i]
    
        orgsum_squares[:,:,i] = orgsum_squares[:,:,i]+ns[i]*(means[:,i]*means[:,i]');
        sum_squares[:,:,i] = sum_squares[:,:,i]+ns[i]*(means[:,i]*means[:,i]');
    
        n = size(baseline[:,labels.==is[i]],2)
        k_n = k_0+n
        v_n = v_0+n
    
        zm_Y = means[:,i]-mu_0
        SS = sum_squares[:,:,i]-n*(means[:,i]*means[:,i]')
        lambda_n = lambda_0 + SS +  k_0*n/(k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-D+1))
    
        log_det_cov[i] = log(det(Sigma))
        orglog_det_cov[i] =  log_det_cov[i]
        inv_cov[:,:,i] = inv(Sigma)
        orginv_cov[:,:,i] =  inv_cov[:,:,i]
    end
    
    dists=Array(Float64,sources)
    for i=1:N
        
        y = datas[:,i]
    
        for j=1:length(uniqs)
            dists[j] =  sqrt(sum((y-orgmeans[:,j]).^2));
        end
    
        choose=find(rand().>=(1-cumsum(dists./sum(dists))))[1]
        allcounts[choose]=allcounts[choose]+1
        counts[choose]=counts[choose]+1
        class_id[i]=choose;
        
        n = allcounts[choose]
        yyT[:,:,i] = y*y';
        
        means[:,choose] = means[:,choose] + (1/n)*(y-means[:,choose]);
        sum_squares[:,:,choose] = sum_squares[:,:,choose] + yyT[:,:,i];
        
        
        k_n = k_0+n;
        v_n = v_0+n;
        
        zm_Y = means[:,choose]-mu_0
        SS = sum_squares[:,:,choose]-n*(means[:,choose]*means[:,choose]')
        lambda_n = lambda_0 + SS + k_0*n/(k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-D+1))
    
        log_det_cov[choose] = log(det(Sigma))
        inv_cov[:,:,choose] = inv(Sigma)
    
    end

    K_plus = max(class_id);
d2 = D/2;
    
    tic()
    totaltime=0
    
    ## start MCMC ---- 
    for iter=1:num_iters
        
        # calculate P for each individual under prior alone
        p_under_prior_alone = p_for_1(pc_gammaln_by_2,pc_log,pc_log_pi,lambda_0,mu_0,k_0,v_0,D,N,datas,p_under_prior_alone)
        
        # run split-merge bit
        
        if iter>20
            (class_id[:,iter],K_plus,sum_squares,means,inv_cov,log_det_cov,counts) = split_merge(datas,class_id[:,iter-1],pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,N,k_0,mu_0,v_0,lambda_0,D,sum_squares,yyT,means,inv_cov,log_det_cov,counts,K_plus,alpha,p_under_prior_alone)
        end 
   
        # run gibbs bit
        
        if iter!=1  && iter.<=20
            class_id[:,iter] = class_id[:,iter-1];
        end

        (class_id[:,iter],K_plus,sum_squares,means,inv_cov,log_det_cov,counts) = crp_gibbs(datas,iter,class_id[:,iter],pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,N,k_0,mu_0,v_0,lambda_0,D,sum_squares,yyT,means,inv_cov,log_det_cov,counts,K_plus,alpha,p_under_prior_alone)

        # update alpha

        alpha = update_alpha(alpha,N,K_plus,a_0,b_0)

        # update prior
        
        (k_0,mu_0) = update_prior(D,K_plus,means,sum_squares,counts,mu_0,v_0,k_0,lambda_0)
        
        # save parameter values
        
        if mod(iter,thin)==0
            K_record[iter/thin] = K_plus
            alpha_record[iter/thin] = alpha
            k_0s[iter/thin]=k_0
            class_ids[:,int(iter/thin)]=class_id[:,iter]
            
        end

  
        # timer

        time_1_iter = toq();
        totaltime = disptime(totaltime,time_1_iter,iter,thin,num_iters,K_record)
        tic()
        
    end


    return(class_ids,k_0s,K_record,alpha_record)

end
