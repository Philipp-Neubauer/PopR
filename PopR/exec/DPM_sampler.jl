function  DPM_sampler(inp::Array{Any,1})

datas = inp[1]
num_iters = inp[2]
thin = inp[3]
Stats = inp[4]
priors = inp[5]
consts = inp[6]
    
    #const (D, N) = size(datas)
    #const pc_max_ind=int(1e4)

    
    #load("define_types.jl")
    require("Julia_code_support.jl")
    require("gibbs_crp.jl")
    require("split-merge.jl")
      
   
    
    alpha = 10

    # initialize structures

    max_class_id = N

    nit=int(num_iters/thin)

    k_0s=zeros(Float,nit)
    class_ids = zeros(Int,(N,nit))
    K_record = zeros(Int,nit)
    alpha_record = zeros(Float,nit)
      
    class_id = zeros(Int,N)
    class_id[1,1] = 1
    counts = zeros(Int,max_class_id,1)
    counts[1] = 1
    K_plus = 1
    p_under_prior_alone = zeros(Float,N);
    
    
    yyT = zeros(Float,(D,D,N))
    for i=1:N
        yyT[:,:,i]=datas[:,i]*datas[:,i]'
    end

    Stats.means[:,1] = datas[:,1]
    Stats.sum_squares[:,:,1] = yyT[:,:,1]
    
    Stats[1] = getlik(consts,priors,Stats[1],datas[:,1],1,false,true)
     
    tic()
    totaltime=0
    
    ## start MCMC ---- 
    for iter=1:num_iters
        
        # calculate P for each individual under prior alone
        p_under_prior_alone = p_for_1(consts,priors,N,datas,p_under_prior_alone)
        
    # if iter==1 || mod(iter,10)==0
        (class_id,K_plus,Stats,counts) = crp_gibbs(datas,iter,class_id,consts,N,priors,yyT,Stats,counts,K_plus,alpha,p_under_prior_alone)
    # end
        
          # run split-merge bit
     #   if iter<(num_iters/10)
        (class_id,K_plus,Stats,counts) = split_merge(datas,class_id,consts,N,priors,yyT,Stats,counts,K_plus,alpha,p_under_prior_alone)
       # end
         
        # update alpha

        alpha = update_alpha(alpha,N,K_plus,priors.a_0,priors.b_0)

        # update prior

        (priors.k_0,priors.mu_0) = update_prior(consts,K_plus,Stats,counts,priors)
        
        # save parameter values
        
        if mod(iter,thin)==0
            K_record[iter/thin] = K_plus
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
