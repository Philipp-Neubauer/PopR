function split_merge(data,class_idd,consts,N,priors,yyT,Stats,counts,K_plus,alpha,p_under_prior_alone)

        classids_temp = deepcopy(class_idd)
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
           
            sstats = NORM(zeros(Float64,D,2),zeros(Float64,D,D,2),zeros(Float64,D,D,2),zeros(Float64,2))
                      
            cstats = NORM( zeros(Float64,D,1),zeros(Float64,D,D),zeros(Float64,D,D),zeros(Float64,1))

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

                cstats,clikelihood = getlik(consts,priors,cstats,y_k,cns,true,true)
                
                clik=clik+clikelihood
                cns=cns+1

                cstats = update_Stats(cstats,y_k,cns,1)
                             
                # calculate individual set likelihoods
                
                if o.==1 || o.==2
                    
                    setlik=setlik+p_under_prior_alone[ind[o]]
                    sstats[o] = update_Stats(sstats[o],y_k,n_S[o],1)
                    continue
                end
  
                for ell = 1:2

                    sstats[ell],likelihood[ell] = getlik(consts,priors,sstats[ell],y_k,n_S[ell],true,true)
                                    
                end
                
                likelihoods = exp(likelihood)

                prob_i[1] = (n_S[1]*likelihoods[1])/sum(n_S'.*likelihoods) # the proba of choosing S_i for individual k
                prob_i[2] = 1-prob_i[1]
                
                setlik=setlik+likelihood[classhelp[o]]
                cprod=cprod*prob_i[classhelp[o]]
                
                n_S[classhelp[o]] = n_S[classhelp[o]]+1

                sstats[classhelp[o]] = update_Stats(sstats[classhelp[o]],y_k,n_S[classhelp[o]],1)
                                              
            end
            
            M_H_prior = exp(lgamma(n_j+n_i)-(lgamma(n_S[1])+lgamma(n_S[2])))/alpha
            M_H_Lik =exp(clik-setlik)
            M_H_rat = M_H_prior*(M_H_Lik)*cprod
            
            if rand().<M_H_rat # accept ?
                #println("accept merge")
                # first update suff-stats of new merged group
                
                counts[c_j] = n_i + n_j
               
                Stats[c_j] = getlik(consts,priors,cstats,0,counts[c_j],false,true)
        
                # then delete old table
                
                class_idd=classids_temp
                K_plus = K_plus-1
                
                hits = [1:c_i-1, c_i+1:(K_plus+1)]
                
                counts[1:K_plus] = counts[hits]
                counts[K_plus+1] = 0
                
                Stats[1:K_plus] =  Stats[hits]
                Stats[K_plus+1] =  0
                
            end
            
        else # split - this is essentially the same thing only that the M-H is slightly different
            
            # reallocate, keeping track of the probas for the M-H ratio
            # for this I need the sufficient stats for the 'new' groups in
            # the merge case starting from one single obs in each group
            # [RANDOM SAMPLE THE ORDER! [Dahl 2003]]
            n_S=ones(Int16,2,1) # keeps track of the number of items at each 'new-old' component
           
            sstats = NORM(zeros(Float64,D,2),zeros(Float64,D,D,2),zeros(Float64,D,D,2),zeros(Float64,2))
                      
            cstats = NORM( zeros(Float64,D,1),zeros(Float64,D,D),zeros(Float64,D,D),zeros(Float64,1))
            
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
                
                cstats,clikelihood = getlik(consts,priors,cstats,y_k,cns,true,true)
                
                clik=clik+(clikelihood)
                cns=cns+1

                cstats = update_Stats(cstats,y_k,cns,1)
                
                # calculate individual set likelihoods
                
                if o.==1 || o.==2
                    
                    setlik=setlik+p_under_prior_alone[ind[o]]
                    sstats[o] = update_Stats(sstats[o],y_k,n_S[o],1)
                    continue
                end
                
                for ell = 1:2

                    sstats[ell],likelihood[ell] = getlik(consts,priors,sstats[ell],y_k,n_S[ell],true,true)
                                    
                end
                
                likelihoods = exp(likelihood)
                prob_i[1] = (n_S[1]*likelihoods[1])/sum(n_S'.*likelihoods) # the proba of choosing S_i for individual k
                prob_i[2] = 1-prob_i[1]
                
                if rand().<prob_i[1] # S_i is chosen, this time I need to keep track of which one k was allocated to                    
                    setlik=setlik+likelihood[1]
                    cprod=cprod*prob_i[1]
                    
                    classids_temp[k]=K_plus+1
                    n_S[1] = n_S[1]+1
                    
                    sstats[1] = update_Stats(sstats[1],y_k,n_S[1],1)
                    
                else # S_j is chosen
                    
                    setlik=setlik+likelihood[2]
                    cprod=cprod*prob_i[2]
                    
                    n_S[2] = n_S[2]+1
                    
                    sstats[2] = update_Stats(sstats[2],y_k,n_S[2],1)
                    
                end
                
            end
            
            M_H_prior = exp((lgamma(n_S[1])+lgamma(n_S[2]))-lgamma(n_j))*alpha
            M_H_Lik = exp(setlik-clik)
            M_H_rat = M_H_prior*(M_H_Lik)*(1/cprod)
            
            if rand().<M_H_rat #&& any(classids_temp!=class_idd) # accept ?
                # println("accept split")
                # first update suff-stats of new groups
                counts[c_j] = n_S[2]
                counts[K_plus+1] = n_S[1]

                Stats[c_j] = getlik(consts,priors,sstats[2],0,counts[c_j],false,true)
                Stats[K_plus+1] = getlik(consts,priors,sstats[1],0,counts[K_plus+1],false,true)
                
                K_plus=K_plus+1
                class_idd=copy(classids_temp)
            end
            
        end

       return(class_idd,K_plus,Stats,counts)
        
end
