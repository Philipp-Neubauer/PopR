function  DPM_sampler_fix(num_iters,thin,Stats,priors,consts)

    load("Julia_code_support.jl")
    load("gibbs_crp_fix.jl")
    load("split-merge_fix.jl")

     alpha = 1

    # initialize structures

    nit=int(num_iters/thin)
    k_0s=Array(Float64,nit)
    K_plus = consts.sources
    counts = zeros(Int64,consts.N,1)
    allcounts = zeros(Int64,consts.N,1)
    class_id = Array(Int64,consts.N)
    class_ids = Array(Int64,(consts.N,nit))
    K_record = Array(Int64,nit)
   
    alpha_record = Array(Float64,nit)
    p_under_prior_alone = Array(Float64,consts.N);

    consts,Stats,allcounts = preallocate(consts,Stats,priors,allcounts)   
    Stats,allcounts,counts = preclass(consts,Stats,priors,allcounts,counts,class_id)
    tic()
    totaltime=0
    
    ## start MCMC ---- 
   for iter=1:num_iters
        
        # calculate P for each individual under prior alone
        p_under_prior_alone = p_for_1(consts,priors,consts.N,consts.datas,p_under_prior_alone)
        
        # run split-merge bit
        
        if iter>10
            (class_id,K_plus,Stats,counts,allcounts) = split_merge(class_id,consts,priors,Stats,counts,allcounts,K_plus,alpha,p_under_prior_alone)
         end
         
        # run gibbs bit      
      
        (class_id,K_plus,Stats,counts,allcounts) = crp_gibbs(class_id,consts,priors,Stats,allcounts,counts,K_plus,alpha,p_under_prior_alone)

        # assert(all(round((Stats.means)*counts/90,4) .==round( mean(datas,2),4)))

        # update alpha

        alpha = update_alpha(alpha,N,sum(counts.!=0),priors.a_0,priors.b_0)

        # update prior
        
        (priors.k_0,priors.mu_0) = update_prior(consts,K_plus,Stats,counts,priors)
        
        # save parameter values
        
        if mod(iter,thin)==0
            K_record[iter/thin] = sum(counts.!=0)
            alpha_record[iter/thin] = alpha
            k_0s[iter/thin]=priors.k_0
            class_ids[:,int(iter/thin)]=class_id
            
        end
  
        # timer

        time_1_iter = toq();
        totaltime = disptime(totaltime,time_1_iter,iter,thin,num_iters,K_record)
        tic()
        
    end


    return(class_ids,k_0s,K_record,alpha_record)

end
