
    
        type STUD
            
            pc_max_ind::Int
            pc_gammaln_by_2::Array{Float,1}
            pc_log_pi::Float
            pc_log::Array{Float,1}
            D::Int
        end

        # precompute student-T posterior predictive distribution constants
        consts = STUD(pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D)
        
   
  
  

   
    # prior composite type
    type MNIW
        
        a_0
        b_0
        k_0
        v_0
        mu_0
        lambda_0
        
    end
    
    priors= MNIW(single_priors[1],single_priors[2],single_priors[3],single_priors[4],single_priors[5:end],matrix_priors)
    
    # stats composite type
    type NORM
        
        means
        sum_squares
        inv_cov
        log_det_cov
        
    end
    
     # set up NORM type stats
    Stats=NORM(zeros(Float,(D,N)),zeros(Float,(D,D,N)),Array(Float,(D,D,N)),Array(Float,(N)))
    
    # define how to access subsets of individuals
    
    function ref(A::NORM,k::Any)

        NORM(A.means[:,k],A.sum_squares[:,:,k],A.inv_cov[:,:,k],A.log_det_cov[k])
            
    end

     # define how to assign subsets of individuals
    
    function assign(A::NORM,B::NORM,k::Any)

        A.means[:,k]=B.means
        A.sum_squares[:,:,k] = B.sum_squares
        A.inv_cov[:,:,k] = B.inv_cov
        A.log_det_cov[k] = B.log_det_cov
    
    end

    function assign(A::NORM,B::Number,k::Any)

        A.means[:,k]=B
        A.sum_squares[:,:,k] = B
        A.inv_cov[:,:,k] = B
        A.log_det_cov[k] = B
    
    end
    
#end


