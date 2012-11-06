
    for i=1:N
    
        if iter==1 && i==1
            continue
        end
       
        y = datas[:,i]
       
        old_class_ids= class_idd[i]
        tempStats = deepcopy(Stats);

        if old_class_ids != 0
            counts[old_class_ids] = counts[old_class_ids] -1
                       
            if counts[old_class_ids].==0
                # delete class compact all data structures
                
                hits = class_idd.>=old_class_ids;
                class_idd[hits] = class_idd[hits]-1
                K_plus = K_plus-1
                
                hits = [1:old_class_ids-1, old_class_ids+1:(K_plus+1)]
                counts[1:K_plus] = counts[hits];
                counts[K_plus+1] = 0;           

                
                tempStats[1:K_plus] =  tempStats[hits]
                tempStats[K_plus+1] =  0
                                             
            else
                
                 tempStats[old_class_ids]=update_Stats(tempStats[old_class_ids],y,counts[old_class_ids],-1)

            end
        end
        
        # prior with new source prob.
        if iter != 1
            prior = [counts[1:K_plus]; alpha]/(N-1+alpha)
        else
            prior = [counts[1:K_plus]; alpha]/(i-1+alpha)
        end
        
        likelihood = zeros(Float64,(length(prior),1))
              
        for ell = 1:K_plus
           
            if old_class_ids.== ell || old_class_ids.== 0
                (tempStats[ell],likelihood[ell]) = getlik(consts,priors,tempStats[ell],y,counts[ell],true,true)
                
            else
                likelihood[ell] = getlik(consts,priors,tempStats[ell],y,counts[ell],true,false)
            end
           
            
        
        end
        
        likelihood[K_plus+1] = p_under_prior_alone[i]
        
        likelihood = exp(likelihood-max(likelihood))

        likelihood = likelihood/sum(likelihood);
        
        # compute the posterior over seating assignment for datum i
        posterior = prior.*likelihood; # this is actually a proportionality
        # normalize the posterior
        posterior = posterior/sum(posterior);
        
        # pick the new source
        cdf = cumsum(posterior);
        rn = rand()
        
        new_class_ids = find(cdf.>rn)[1]
        
        counts[new_class_ids] = counts[new_class_ids]+1
               
        if new_class_ids .== (K_plus+1)
           K_plus = K_plus+1;
        end
        
        if old_class_ids .!= new_class_ids
          
           tempStats[new_class_ids]=update_Stats(tempStats[new_class_ids],y,counts[new_class_ids],1)

           tempStats[new_class_ids] = getlik(consts,priors,tempStats[new_class_ids],y,counts[new_class_ids],false,true)
           Stats = deepcopy(tempStats) 
        end
        # record the new source
        class_idd[i] = new_class_ids

assert(all(round(mean(datas[:,class_idd.==class_idd[i]],2),4).==round(Stats[class_idd[i]].means,4)))
        
    end
