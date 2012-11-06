function split_merge(data,class_idd,pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,N,k_0,mu_0,v_0,lambda_0,D,sum_squares,yyT,means,inv_cov,log_det_cov,counts,K_plus,alpha,p_under_prior_alone)

        classids_temp = copy(class_idd)
        # choose individuals at random
        ind = randi(N,2)
        c_i = classids_temp[ind[1]]
        c_j = classids_temp[ind[2]]
        
        n_i = sum(classids_temp.==c_i)
        n_j = sum(classids_temp.==c_j)
        
        if c_i !=c_j
            ## merge
            classids_temp[classids_temp.==c_i] = c_j# rename class_idd of ind[1] to that of ind[2]
            hits = classids_temp .>= c_i # get class_idds of id_s greater than the merged one
            classids_temp[hits] = classids_temp[hits]-1 # compact classids_temp
            
            # backwards reallocate, keeping track of the probas for the M-H ratio
            # for this I need the sufficient stats for the 'old' groups in
            # the merge case starting from one single obs in each group
            # (RANDOM SAMPLE THE ORDER! [Dahl 2003]]
            n_S=ones(Int16,2,1) # keeps track of the number of items at each 'new-old' component
           
            meens = zeros(Float64,D,2)
            meens[:,1] = datas[:,ind[1]]
            meens[:,2] = datas[:,ind[2]]

            sum_s = zeros(Float64,D,D,2)
            sum_s[:,:,1] = yyT[:,:,ind[1]]
            sum_s[:,:,2] = yyT[:,:,ind[2]]
            
            cmeans = zeros(Float64,D,1)
            csums = zeros(Float64,D,D)
            cns = 0
            clik=0;cprod=1;setlik=0;likelihood=zeros(Float64,(1,2))
            prob_i=zeros(Float64,(2,1))

            idxar = Array(Bool,N,1)
            for k=1:N
                idxar[k] =class_idd[k].==c_j || class_idd[k].==c_i
            end
            idxar[ind]=0
            ixxs = find(idxar)
            ixxs=[ind[1]; ind[2]; ixxs]

            classhelp = ones(Int16,length(ixxs))
            classhelp[class_idd[ixxs].==c_i]=2
            
            for o = [1, 2, randperm(n_i+n_j-2)+2]
                
                k=ixxs[o]
                
                # calculate the combined likelihood
                
                y_k = datas[:,k]

                clikelihood = student_lp(pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,y_k,cns,cmeans,k_0,mu_0,v_0,lambda_0,D,csums,0)
                
                clik=clik+clikelihood
                cns=cns+1
                cmeans = cmeans + (1/cns)*(y_k-cmeans)
                csums = csums + yyT[:,:,k]
                
                # calculate individual set likelihoods
                
                if o.==1 || o.==2
                    
                    setlik=setlik+p_under_prior_alone[ind[o]]
                    
                    continue
                end
  
                for ell = 1:2

                    likelihood[ell] = student_lp(pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,y_k,n_S[ell],meens[:,ell],k_0,mu_0,v_0,lambda_0,D,sum_s[:,:,ell],0)[1]
                                    
                end
                
                likelihoods = exp(likelihood)

                prob_i[1] = (n_S[1]*likelihoods[1])/sum(n_S'.*likelihoods) # the proba of choosing S_i for individual k
                prob_i[2] = 1-prob_i[1]
                
                setlik=setlik+likelihood[classhelp[o]]
                cprod=cprod*prob_i[classhelp[o]]
                
                n_S[classhelp[o]] = n_S[classhelp[o]]+1
                
                meens[:,classhelp[o]] = meens[:,classhelp[o]]+ (1/n_S[classhelp[o]])*(y_k-meens[:,classhelp[o]])
                sum_s[:,:,classhelp[o]] = sum_s[:,:,classhelp[o]] + yyT[:,:,k]
                               
            end
            
            M_H_prior = exp(lgamma(n_j+n_i)-(lgamma(n_S[1])+lgamma(n_S[2])))/alpha
            M_H_Lik =exp(clik-setlik)
            M_H_rat = M_H_prior*(M_H_Lik)*cprod
            
            if rand().<M_H_rat[1] # accept ?
                #println("accept merge")
                # first update suff-stats of new merged group
                
                counts[c_j] = n_i + n_j
                means[:,c_j] = (mean(datas[:,class_idd.==c_j],2)*n_j+n_i*mean(datas[:,class_idd.==c_i],2))/(n_i+n_j)
                sum_squares[:,:,c_j] = sum_squares[:,:,c_j] + sum_squares[:,:,c_i]
                
                # update relevant quantities for student -t for c_j
                
                m_Y = means[:,c_j]
                n=counts[c_j]
                k_n = k_0+n
                v_n = v_0+n
                
                S = (sum_squares[:,:,c_j] - n*(m_Y*m_Y'))
                zm_Y = m_Y-mu_0
                lambda_n = lambda_0 + S  + k_0*n/(k_0+n)*zm_Y*zm_Y'
                Sigma = (lambda_n*(k_n+1)/(k_n*(v_n-D+1)))'
                
                log_det_cov[c_j] = log(det(Sigma))
                inv_cov[:,:,c_j] = Sigma^-1
                
                # then delete old table
                
                class_idd=classids_temp
                K_plus = K_plus-1
                
                hits = [1:c_i-1, c_i+1:(K_plus+1)]
                means[:,1:K_plus] = means[:,hits]
                means[:,K_plus+1] = 0
                sum_squares[:,:,1:K_plus] = sum_squares[:,:,hits]
                sum_squares[:,:,1+K_plus] = 0
                counts[1:K_plus] = counts[hits]
                counts[K_plus+1] = 0
                
                log_det_cov[1:K_plus] = log_det_cov[hits]
                log_det_cov[K_plus+1] = 0
                inv_cov[:,:,1:K_plus] = inv_cov[:,:,hits]
                inv_cov[:,:,K_plus+1] = 0
                
            end
            
        else # split - this is essentially the same thing only that the M-H is slightly different
            
            # reallocate, keeping track of the probas for the M-H ratio
            # for this I need the sufficient stats for the 'new' groups in
            # the merge case starting from one single obs in each group
            # [RANDOM SAMPLE THE ORDER! [Dahl 2003]]
            n_S=ones(Int16,2,1) # keeps track of the number of items at each 'new-old' component
           
            meens = zeros(Float64,D,2)
            meens[:,1] = datas[:,ind[1]]
            meens[:,2] = datas[:,ind[2]]

            sum_s = zeros(Float64,D,D,2)
            sum_s[:,:,1] = yyT[:,:,ind[1]]
            sum_s[:,:,2] = yyT[:,:,ind[2]]
            
            cmeans = zeros(Float64,(D,1))
            csums = zeros(Float64,(D,D))
            cns = 0
            clik=0;cprod=1;setlik=0;likelihood=zeros(Float64,(1,2))
            prob_i=zeros(Float64,(2,1))
            
            classids_temp[ind[1]]=K_plus+1 # ind 1 marks new class
                 
            ixxs = find(class_idd.==c_j)
            ixxs=ixxs[ixxs.!=ind[1]]
            ixxs=ixxs[ixxs.!=ind[2]]
            ixxs=[ind[1], ind[2], ixxs]
            
            for o = [1, 2, randperm(n_j-2)+2]
                
                k=ixxs[o]
                
                # calculate the combined likelihood
                y_k = datas[:,k]
                
                clikelihood = student_lp(pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,y_k,cns,cmeans,k_0,mu_0,v_0,lambda_0,D,csums,0)
                
                clik=clik+(clikelihood)
                cns=cns+1
                cmeans = cmeans + (1/cns)*(y_k-cmeans)
                csums = csums + yyT[:,:,k]
                
                # calculate individual set likelihoods
                
                if o.==1 || o.==2
                    
                    setlik=setlik+p_under_prior_alone[ind[o]]
                    
                    continue
                end
                
                for ell = 1:2

                    likelihood[ell] = student_lp(pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,y_k,n_S[ell],meens[:,ell],k_0,mu_0,v_0,lambda_0,D,sum_s[:,:,ell],0)[1]
                                    
                end
                
                likelihoods = exp(likelihood)
                prob_i[1] = (n_S[1]*likelihoods[1])/sum(n_S'.*likelihoods) # the proba of choosing S_i for individual k
                prob_i[2] = 1-prob_i[1]
                
                if rand().<prob_i[1] # S_i is chosen, this time I need to keep track of which one k was allocated to
                    
                    setlik=setlik+likelihood[1]
                    cprod=cprod*prob_i[1]
                    
                    classids_temp[k]=K_plus+1
                    n_S[1] = n_S[1]+1
                    
                    meens[:,1] = meens[:,1]+ (1/n_S[1])*(y_k-meens[:,1])
                    sum_s[:,:,1] = sum_s[:,:,1] + yyT[:,:,k]
                    
                else # S_j is chosen
                    
                    setlik=setlik+likelihood[2]
                    cprod=cprod*prob_i[2]
                    
                    n_S[2] = n_S[2]+1
                    
                    meens[:,2] = meens[:,2]+ (1/n_S[2])*(y_k-meens[:,2])
                    sum_s[:,:,2] = sum_s[:,:,2] + yyT[:,:,k]
                    
                end
                
            end
            
            M_H_prior = exp((lgamma(n_S[1])+lgamma(n_S[2]))-lgamma(n_j))*alpha
            M_H_Lik = exp(setlik-clik)
            M_H_rat = M_H_prior*(M_H_Lik)*(1/cprod)
            
            if rand().<M_H_rat[1] #&& any(classids_temp!=class_idd) # accept ?
                # println("accept split")
                # first update suff-stats of new groups
                counts[c_j] = n_S[2]
                counts[K_plus+1] = n_S[1]
                means[:,c_j] = mean(datas[:,classids_temp.==c_j],2)
                means[:,K_plus+1] = mean(datas[:,classids_temp.==(K_plus+1)],2)
                
                sum_squares[:,:,K_plus+1] = reshape(sum(yyT[:,:,classids_temp.==(K_plus+1)],3),D,D)
                sum_squares[:,:,c_j] = sum_squares[:,:,c_j] - sum_squares[:,:,K_plus+1]
                
                # update relevant quantities for student -t for c_j and c_i
                
                m_Y = means[:,c_j]
                n=counts[c_j]
                k_n = k_0+n
                v_n = v_0+n
                
                S = (sum_squares[:,:,c_j] - n*(m_Y*m_Y'))
                zm_Y = m_Y-mu_0
                lambda_n = lambda_0 + S  +  k_0*n/(k_0+n)*(zm_Y*zm_Y')
                Sigma = (lambda_n*(k_n+1)/(k_n*(v_n-D+1)))'
                
                log_det_cov[c_j] = log(det(Sigma))
                inv_cov[:,:,c_j] = Sigma^-1
                
                m_Y = means[:,K_plus+1]
                n=counts[K_plus+1]
                k_n = k_0+n
                v_n = v_0+n
                
                S = (sum_squares[:,:,K_plus+1] - n*(m_Y*m_Y'))
                zm_Y = m_Y-mu_0
                lambda_n = lambda_0 + S  +  k_0*n/(k_0+n)*zm_Y*zm_Y'
                Sigma = (lambda_n*(k_n+1)/(k_n*(v_n-D+1)))'
                
                log_det_cov[K_plus+1] = log(det(Sigma))
                inv_cov[:,:,K_plus+1] = Sigma^-1
                
                K_plus=K_plus+1
                class_idd=copy(classids_temp)
            end
            
        end

       return(class_idd,K_plus,sum_squares,means,inv_cov,log_det_cov,counts)
        
end
