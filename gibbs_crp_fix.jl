function crp_gibbs(datas,iter,class_idd,pc_max_ind,pc_gammaln_by_2,pc_log_pi,pc_log,N,k_0,mu_0,v_0,lambda_0,D,sum_squares,yyT,means,inv_cov,log_det_cov,counts,K_plus,alpha, p_under_prior_alone)

    # take this out later !!!!
    #class_ids=class_id[:,iter]

    for i=1:N
    
        y = datas[:,i]
       
        old_class_ids= class_idd[i]
       
      
        counts[old_class_ids] = counts[old_class_ids] -1
        allcounts[old_class_ids] = allcounts[old_class_ids] -1

        old_class_log_det_Sigma = log_det_cov[old_class_ids]
        old_class_inv_Sigma = inv_cov[:,:,old_class_ids]
            
        if counts[old_class_ids].==0 && old_class_id>sources
                # delete class compact all data structures
                
            hits = class_idd.>=old_class_ids;
            class_idd[hits] = class_idd[hits]-1;
            K_plus = K_plus-1;
                
            hits = [1:(old_class_ids-1), (old_class_ids+1):(K_plus+1)]

            means[:,1:K_plus] = means[:,hits];
            means[:,K_plus+1] = 0;
            sum_squares[:,:,1:K_plus] = sum_squares[:,:,hits];
            sum_squares[:,:,1+K_plus] = 0;
            counts[1:K_plus] = counts[hits];
            counts[K_plus+1] = 0;
            allcounts[1:K_plus] = allcounts[hits];
            allcounts[K_plus+1] = 0;
                
            log_det_cov[1:K_plus] = log_det_cov[hits];
            log_det_cov[K_plus+1] = 0;
            inv_cov[:,:,1:K_plus] = inv_cov[:,:,hits];
            inv_cov[:,:,K_plus+1] = 0;
                
                
        else
            means[:,old_class_ids] = (1/(allcounts[old_class_ids]))*((allcounts[old_class_ids]+1)*means[:,old_class_ids] - y);
            sum_squares[:,:,old_class_ids] = sum_squares[:,:,old_class_ids] - yyT[:,:,i];
        end
        

        priornum = Array(Float64,K_plus)
        # complete the CRP prior with new source prob.
       
        priornum[counts[1:K_plus].!=0] = counts[counts[1:K_plus].!=0]
        priornum[counts[1:K_plus].==0] = repmat(alpha,sum[counts[1:K_plus].==0),1)/(sum(counts[1:K_plus].==0)+1)
        prior = [priornum';alpha/(sum(counts[1:K_plus)==0]+1)]/(N-1+alpha);
               
        likelihood = zeros(Float64,(length(prior),1))
              
        # as per Radford's Alg. 3 compute the posterior predictive
        # probabilities in two scenerios, 1) we will evaluate the
        # likelihood of sitting at all of the existing sources by computing
        # the probability of the datapoint under the posterior predictive
        # distribution with all points sitting at that source considered and
        # 2) we will compute the likelihood of the point under the
        # posterior predictive distribution with no observations
        
        for ell = 1:K_plus
            # get the class ids of the points sitting at source l
            
            n = allcounts[ell]
            
            m_Y = means[:,ell]
            mu_n = k_0/(k_0+n)*mu_0 + n/(k_0+n)*m_Y
            k_n = k_0+n
            v_n = v_0+n
            
            # set up variables for Gelman's formulation of the Student T
            # distribution
            v = v_n-D+1
            mu = mu_n
            
            
            # if old_class_ids == ell means that this point used to sit at
            # source ell, all of the sufficient statistics have been updated
            # in sum_squares, counts, and means but that means that we have
            # to recompute log_det_Sigma and inv_Sigma.  if we reseat the
            # particle at its old source then we can put the old
            # log_det_Sigma and inv_Sigma back, otherwise we need to update
            # both the old source and the new source
           
            
                if old_class_ids .== ell
                    S = (sum_squares[:,:,ell] - n*(m_Y*m_Y'))
                    zm_Y = m_Y-mu_0
                    lambda_n = lambda_0 + S  + k_0*n/(k_0+n)*(zm_Y)*(zm_Y)'
                    Sigma = (lambda_n*(k_n+1)/(k_n*(v_n-D+1)))'
                                       
                    log_det_Sigma = log(det(Sigma))
                    inv_Sigma = (Sigma)^-1
                    log_det_cov[old_class_ids] = log_det_Sigma
                    inv_cov[:,:,old_class_ids] = inv_Sigma
                else
                    log_det_Sigma = log_det_cov[ell]
                    inv_Sigma = inv_cov[:,:,ell];
                end
                       
            vd = v+D;
            d2 = D/2
            # the log likelihood for class ell
            likelihood[ell] = (pc_gammaln_by_2[vd] - (pc_gammaln_by_2[v] + d2*pc_log[v] + d2*pc_log_pi) - .5*log_det_Sigma- (vd/2)*log(1+(1/v)*(y-mu)'*inv_Sigma*(y-mu)))[1]
            
        end
        
        likelihood[K_plus+1] = p_under_prior_alone[i];
        
        likelihood = exp(likelihood-max(likelihood));
        likelihood = likelihood/sum(likelihood);
        
        # compute the posterior over seating assignment for datum i
        posterior = prior.*likelihood; # this is actually a proportionality
        # normalize the posterior
        posterior = posterior/sum(posterior);
        
        # pick the new source
        cdf = cumsum(posterior);
        rn = rand()
        
        new_class_ids = find(cdf.>rn)[1]
        
        counts[new_class_ids] = counts[new_class_ids]+1;
        means[:,new_class_ids] = means[:,new_class_ids]+ (1/counts[new_class_ids])*(y-means[:,new_class_ids]);
        sum_squares[:,:,new_class_ids] = sum_squares[:,:,new_class_ids] + yyT[:,:,i];
        
        if new_class_ids .== (K_plus+1)
           K_plus = K_plus+1;
        end
        
        if old_class_ids .== new_class_ids
            # we don't need to compute anything new as the point was
            # already sitting at that source and the matrix inverse won't
            # change
            log_det_cov[old_class_ids] = old_class_log_det_Sigma;
            inv_cov[:,:,old_class_ids] = old_class_inv_Sigma;
        else
            # the point changed sources which means that the matrix inverse
            # sitting in the old_class_ids slot is appropriate but that the
            # new source matrix inverse needs to be updated
            n = counts[new_class_ids];
            #             if n~=0
            m_Y = means[:,new_class_ids];
            k_n = k_0+n;
            v_n = v_0+n;
            
            # set up variables for Gelman's formulation of the Student T
            # distribution
            S = (sum_squares[:,:,new_class_ids] - n*(m_Y*m_Y'));
            zm_Y = m_Y-mu_0;
            lambda_n = lambda_0 + S  + k_0*n/(k_0+n)*(zm_Y)*(zm_Y)';
            Sigma = (lambda_n*(k_n+1)/(k_n*(v_n-D+1)))';
            
            log_det_cov[new_class_ids] = log(det(Sigma));
            inv_cov[:,:,new_class_ids] = Sigma^-1;
        end
        
        # record the new source
        class_idd[i] = new_class_ids;
        
    end
    return(class_idd,K_plus,sum_squares,means,inv_cov,log_det_cov,counts)
    
    end
