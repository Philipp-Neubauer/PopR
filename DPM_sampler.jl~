function  DPM_sampler(num_iters,thin)

single_priors = dlmread("single_priors.csv",",",Float64)
single_priors=single_priors[2:end,2]

matrix_priors = dlmread("matrix_priors.csv",",",Float64)
matrix_priors=matrix_priors[2:end,2:end]

datas =dlmread("datas.csv",",",Float64)
global datas=datas[2:end,2:end]'

    load("Julia_code_support.jl")
    load("gibbs_crp.jl")
    load("split-merge.jl")

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

    k_0s=zeros(Float64,nit)
    K_plus = 1
    class_id = zeros(Int16,(N,num_iters))
    class_ids = zeros(Int16,(N,nit))
    K_record = zeros(Int16,nit)
    alpha_record = zeros(Float64,nit)
    
    # seat the first customer at the first table
    class_id[1,1] = 1
    
    # precompute student-T posterior predictive distribution constants
    const pc_max_ind = 1e5
    const pc_gammaln_by_2 = lgamma((1:pc_max_ind)/2)
    const pc_log_pi = log(pi)
    const pc_log = log(1:pc_max_ind)
    
    means = zeros(Float64,(D,max_class_id))
    sum_squares = zeros(Float64,(D,D,max_class_id))
    inv_cov = zeros(Float64,(D,D,max_class_id))
    log_det_cov = zeros(Float64,(max_class_id))
    counts = zeros(Int16,max_class_id,1)
    counts[1] = 1
    
    p_under_prior_alone = zeros(Float64,N);
    
    # sit first individual at first table
    y = datas[:,1]
    yyT = y*y'
    
    (ldc,ic) = student_lp(pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,y,1,y,k_0,mu_0,v_0,lambda_0,D,yyT,1)
    
    means[:,1] = y
    sum_squares[:,:,1] = yyT
    counts[1] = 1
    log_det_cov[1] = ldc
    inv_cov[:,:,1] = ic
    
    yyT = zeros(Float64,(D,D,N))
    
    for i=1:N
        yyT[:,:,i]=datas[:,i]*datas[:,i]'
    end
    
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
