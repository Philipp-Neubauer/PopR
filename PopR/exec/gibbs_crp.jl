function crp_gibbs(datas,iter,class_idd,consts,N,priors,yyT,stats,counts,K_plus,alpha,p_under_prior_alone)

    # take this out later !!!!
    #class_idd=class_id
    #class_idd_old=deepcopy(class_idd)
    #stats_old = deepcopy(stats)
    #counts_old=deepcopy(counts)
    #K_plus_old=deepcopy(K_plus)
    
    for i=1:N
    
        if iter==1 && i==1
            continue
        end
       
        y = datas[:,i]
        
       
        old_class_ids= class_idd[i]

        # temporary variables
        tempstats = deepcopy(stats);
        class_idd_temp = class_idd-0
        K_plus_temp = K_plus-0
        tempcounts = counts-0          
        dels=false
        
        if old_class_ids != 0
            
            tempcounts[old_class_ids] = tempcounts[old_class_ids] -1
                       
            if tempcounts[old_class_ids].==0
                dels = true
                # delete class compact all data structures
                               
                hits = class_idd_temp.>=old_class_ids
                class_idd_temp[hits] = class_idd_temp[hits]-1
                K_plus_temp = K_plus-1
                
                hits = [1:old_class_ids-1, old_class_ids+1:(K_plus_temp+1)]
                tempcounts[1:K_plus_temp] = tempcounts[hits];
                tempcounts[K_plus_temp+1] = 0;           
                old_class_ids=K_plus_temp+1
                
                tempstats[1:K_plus_temp] =  tempstats[hits]
                tempstats[K_plus_temp+1] =  0
                                             
            else
                
                tempstats[old_class_ids]=update_stats(tempstats[old_class_ids],y,tempcounts[old_class_ids],-1)
                
            end
        end
        
        # prior with new source prob.
        if iter != 1
            prior = [tempcounts[1:K_plus_temp]; alpha]/(N-1+alpha)
        else
            prior = [tempcounts[1:K_plus_temp]; alpha]/(i-1+alpha)
        end
        
        likelihood = zeros(Float,(length(prior),1))
              
        for ell = 1:K_plus_temp
           
            if old_class_ids.== ell || old_class_ids.== 0
                (tempstats[ell],likelihood[ell]) = getlik(consts,priors,tempstats[ell],y,tempcounts[ell],true,true)
                
            else
                likelihood[ell] = getlik(consts,priors,tempstats[ell],y,tempcounts[ell],true,false)
            end                     
        
        end
        
        likelihood[K_plus_temp+1] = p_under_prior_alone[i]

       # if !all((likelihood.>=0.0) | (likelihood.<=0.0))
        #println(likelihood)
       #     return(class_idd_old,K_plus_old,stats_old,counts_old)
       # end
        
        likelihood = exp(likelihood-max(likelihood))

        likelihood = likelihood/sum(likelihood)
        
        # compute the posterior over seating assignment for datum i
        posterior = prior.*likelihood; # this is actually a proportionality
        # normalize the posterior
        posterior = posterior/sum(posterior);
        
        # pick the new source
        cdf = cumsum(posterior)
        rn = rand()
        
        new_class_ids = find(cdf.>rn)[1]
        
        tempcounts[new_class_ids] = tempcounts[new_class_ids]+1
        newc=false
        if new_class_ids .== (K_plus_temp+1)
        newc=true
           K_plus_temp = K_plus_temp+1;
        end

        # update things either the new class id is != the old class without any re-arrangements, or it's uneuqal to K_plus+1 with re-arragement -- else we can jsut keep everything as is...
        if (old_class_ids .!= new_class_ids && dels==false)||(newc==false && dels==true)

            # record the new source and update variables to values of temporary variables
            class_idd_temp[i] = new_class_ids
            class_idd =class_idd_temp-0
            K_plus = K_plus_temp-0
            counts = tempcounts-0
            
            tempstats[new_class_ids]=update_stats(tempstats[new_class_ids],y,counts[new_class_ids],1)

                    
            tempstats[new_class_ids] = getlik(consts,priors,tempstats[new_class_ids],y,counts[new_class_ids],false,true)
                  
            stats = deepcopy(tempstats) 
        end
    
    end
    
    return(class_idd,K_plus,stats,counts)
    
    end
