#if(Typeof=="N")

    single_priors = dlmread("single_priors.csv",",",Float64)
    single_priors=single_priors[2:end,2]
    
    matrix_priors = dlmread("matrix_priors.csv",",",Float64)
    matrix_priors=matrix_priors[2:end,2:end]
   
    # define normal set types

    pc_max_ind=1e5
    
    if bl.==1

    
        baseline =dlmread("baseline.csv",",",Float64)
        global baseline=baseline[2:end,2:end]'


        label = dlmread("labels.csv",",",Float64)
        label = label[2:end,2]
        label=int(label)

        function unique{T}(A::AbstractArray{T}, sorted::Bool)
            dd = Dict{T, Bool}()
            for a in A dd[a] = true end
            sorted? sort!(keys(dd)): keys(dd)
        end
        
        uniqs=unique(label,true)
        labels=Array(Int64,length(label))
        for i=1:length(uniqs)
            labels[label.==uniqs[i]]=i
        end

        type STUD
            datas
            baseline
            orgsum_squares
            orgmeans
            orginv_cov
            orglog_det_cov
            ns
            labels
            sources
            pc_max_ind
            pc_gammaln_by_2
            pc_log_pi
            pc_log
            D
            N
        end
         
    # precompute student-T posterior predictive distribution constants
    consts = STUD(datas,baseline,Array(Float64,(D,D,max(labels))),Array(Float64,(D,max(labels))),Array(Float64,(D,D,max(labels))),Array(Float64,max(labels)),Array(Int64,max(labels)),labels,max(labels),pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D,N)

    else
    
        type STUD
            
            pc_max_ind
            pc_gammaln_by_2
            pc_log_pi
            pc_log
            D
        end

        # precompute student-T posterior predictive distribution constants
        consts = STUD(pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D)
        
    end
  
  

   
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
    Stats=NORM(zeros(Float64,(D,N)),zeros(Float64,(D,D,N)),Array(Float64,(D,D,N)),Array(Float64,(N)))
    
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


