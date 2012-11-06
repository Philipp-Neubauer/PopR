function create_structs(Typeof,D)

# define structures for normal model
    if(Typeof=="N")

        single_priors = dlmread("single_priors.csv",",",Float64)
        single_priors=single_priors[2:end,2]

        matrix_priors = dlmread("matrix_priors.csv",",",Float64)
        matrix_priors=matrix_priors[2:end,2:end]

       # define normal set types

        type STUD

            pc_max_ind
            pc_gammaln_by_2
            pc_log_pi
            pc_log
            D
        end

        # precompute student-T posterior predictive distribution constants
        consts = STUD(1e5,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D)

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

         Stats=NORM(Array(Float64,(D,max_class_id)),Array(Float64,(D,D,max_class_id)),Array(Float64,(D,D,max_class_id)),Array(Float64,(max_class_id)));

        # define how to access subsets of individuals

        function ref(A::NORM,k::Any)

            NORM(A.means[:,k],A.sum_squares[:,:,k],A.inv_cov[:,:,k],A.log_det_cov[k])
             
         end


    end
    

    # define structures for Multinomial model
    # if(Typeof=="MN")

    
    return consts,priors,Stats

end
