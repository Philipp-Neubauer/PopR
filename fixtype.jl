import Base.getindex
import Base.setindex!
import Base.copy


        baseline=readdlm("baseline.csv",",",Float)[2:end,2:end]'

        label = readdlm("labels.csv",",",Float)[2:end,2]
        label=int(label)  
        

        uniqs=sort(unique(label))
        labels=Array(Int,length(label))
        for i=1:length(uniqs)
            labels[label.==uniqs[i]]=i
        end

        immutable STUD
            datas::Array{Float,2}
            baseline::Array{Float,2}
            orgsum_squares::Array{Float,3}
            orgmeans::Array{Float,2}
            orginv_cov::Array{Float,3}
            orglog_det_cov::Array{Float,1}
            ns::Array{Int,1}
            labels::Array{Int,1}
            sources::Int
            pc_max_ind::Int
            pc_gammaln_by_2::Array{Float,1}
            pc_log_pi::Float
            pc_log::Array{Float,1}
            D::Int
            N::Int
        end
         
    # precompute student-T posterior predictive distribution constants
    consts = STUD(datas,baseline,Array(Float,(D,D,max(labels))),Array(Float,(D,max(labels))),Array(Float,(D,D,max(labels))),Array(Float,max(labels)),Array(Int,max(labels)),labels,int(max(labels)),pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D,N)


     # prior composite type
    type MNIW

        a_0::Float
        b_0::Float
        k_0::Float
        ak_0::Float
        bk_0::Float
        v_0::Float
        mu_0::Array{Float,1}
        lambda_0::Array{Float,2}
        
    end
    
    priors= MNIW(single_priors[1],single_priors[2],single_priors[3],single_priors[4],single_priors[5],single_priors[6],single_priors[7:end],matrix_priors)
    
    # stats composite type
    type NORM
        
        means
        sum_squares
        inv_cov
        log_det_cov
        
    end
    
     # set up NORM type stats
    Stats=NORM(zeros(Float,(D,N)),zeros(Float,(D,D,N)),Array(Float,(D,D,N)),Array(Float,(N)))

# define how to copy
    function copy(A::NORM)

        NORM(A.means,A.sum_squares,A.inv_cov,A.log_det_cov)
            
    end

    # define how to access subsets of individuals
    
    function getindex(A::NORM,k::Int)

        NORM(A.means[:,k],A.sum_squares[:,:,k],A.inv_cov[:,:,k],A.log_det_cov[k])
            
    end

 # define how to access subsets of individuals
    
    function getindex(A::NORM,k::Array{Int,1})

        NORM(A.means[:,k],A.sum_squares[:,:,k],A.inv_cov[:,:,k],A.log_det_cov[k])
            
    end

     # define how to assign subsets of individuals
    
    function setindex!(A::NORM,B::NORM,k::Int)

        A.means[:,k]=B.means
        A.sum_squares[:,:,k] = B.sum_squares
        A.inv_cov[:,:,k] = B.inv_cov
        A.log_det_cov[k] = B.log_det_cov
    
    end

    function setindex!(A::NORM,B::NORM,k::Range1{Int})

        A.means[:,k]=B.means
        A.sum_squares[:,:,k] = B.sum_squares
        A.inv_cov[:,:,k] = B.inv_cov
        A.log_det_cov[k] = B.log_det_cov
    
    end
    
    function setindex!(A::NORM,B::Number,k::Any)

        A.means[:,k]=B
        A.sum_squares[:,:,k] = B
        A.inv_cov[:,:,k] = B
        A.log_det_cov[k] = B
    
    end


