function split_merge(class_idd,consts,priors,Stats,counts,allcounts,K_plus,alpha,p_under_prior_alone)

   classids_temp = copy(class_idd)
    # choose individuals at random
     ind = rand(1:consts.N,2)
    while ind[1]==ind[2]
        ind[2] = rand(1:consts.N,1)[1]
    end
    
    c_i = classids_temp[ind[1]]
    c_j = classids_temp[ind[2]]
        
    n_i = sum(classids_temp.==c_i)
    n_j = sum(classids_temp.==c_j)
    cs=[c_i,c_j]
        
    if c_i !=c_j
        ## merge
           
        if (c_i>consts.sources && c_j>consts.sources) || ((c_i<=consts.sources && c_j<=consts.sources) && n_i >= n_j) || (c_i<=consts.sources && c_j>consts.sources)
                
            classids_temp[class_idd.==c_j] = c_i# rename class_id of ind(1) to that of ind(2)
                
            if c_j>consts.sources || c_i>consts.sources # only need to compact if one is not from baseline
                    
                hits = classids_temp .>= c_j; # get class_ids of id_s greater than the merged one
                classids_temp[hits] = classids_temp[hits]-1 # compact class_id_temp
                    
            end
                
        else
                
            classids_temp[class_idd.==c_i] = c_j;# rename class_id of ind(1) to that of ind(2)
                
            if c_j>consts.sources || c_i>consts.sources # only need to compact if one is not from baseline
                    
                hits = classids_temp .>= c_i; # get class_ids of id_s greater than the merged one
                classids_temp[hits] = classids_temp[hits]-1  # compact class_id_temp
                
            end
                
        end
                   
            
            # backwards reallocate, keeping track of the probas for the M-H ratio
            # for this I need the sufficient stats for the 'old' groups in
            # the merge case starting from one single obs in each group
            # (RANDOM SAMPLE THE ORDER! [Dahl 2003]]
        n_S=ones(Int16,2,1) # keeps track of the number of items at each 'new-old' component
           
        sstats = NORM(zeros(Float,D,2),zeros(Float,D,D,2),zeros(Float,D,D,2),zeros(Float,2))
                      
        cstats = NORM( zeros(Float,D),zeros(Float,D,D),zeros(Float,D,D),zeros(Float))

        if c_i<=consts.sources || c_j<=consts.sources

                
            inds =cs[cs.<=consts.sources]

            a=0
            for i=inds
                a=a+1
                n_S[a] = n_S[a]+consts.ns[i]  
                sstats.means[:,a] = consts.orgmeans[:,i]
                sstats.sum_squares[:,:,a] = consts.orgsum_squares[:,:,i]
            end
                

        end

        if (c_i>consts.sources && c_j>consts.sources)
            cns = 0  
            merga = c_i;
            cnt=[1, 2]; # which way around do we start ?
            antimerga = c_j;
            nmerge = n_j;                
                
        elseif  (c_i<=consts.sources && c_j>consts.sources)
            cns = consts.ns[c_i]
            merga = c_i;
            cnt=[1, 2]; # which way around do we start ?                
            antimerga = c_j;
            nmerge = n_j;
            cstats.means = consts.orgmeans[:,c_i]
            cstats.sum_squares = consts.orgsum_squares[:,:,c_i]

        elseif  (c_j<=consts.sources && c_i>consts.sources)
            cns = consts.ns[c_j] 
            merga = c_j;
            cnt=[2,1];
            antimerga = c_i;
            nmerge = n_i;
            cstats.means = consts.orgmeans[:,c_j]
            cstats.sum_squares = consts.orgsum_squares[:,:,c_j]
            
        elseif   (c_i<=consts.sources && c_j<=consts.sources)
       
            if n_i >= n_j
                merga = c_i;
                cnt=[1, 2]; # which way around do we start ?                
                antimerga = c_j;
                nmerge = n_j;
            else 
                merga = c_j;
                cnt=[2,1];
                antimerga = c_i;
                nmerge = n_i;
            end
            cns = consts.ns[c_i]+consts.ns[c_j]
            cstats.means = (consts.orgmeans[:,c_j]*consts.ns[c_j]+consts.orgmeans[:,c_i]*consts.ns[c_i])/cns
            cstats.sum_squares = consts.orgsum_squares[:,:,c_j]+consts.orgsum_squares[:,:,c_i]
        
        end
            
        
        clik=0;cprod=1;setlik=0;likelihood=zeros(Float,(1,2))
        prob_i=zeros(Float,(2,1))

        idxar = Array(Bool,N,1)
        for k=1:N
            idxar[k] =class_idd[k].==c_j || class_idd[k].==c_i
        end
        idxar[ind]=0
        ixxs = find(idxar)
        ixxs=[ind[1]; ind[2]; ixxs]

        classhelp = ones(Int,length(ixxs))
        classhelp[class_idd[ixxs].==c_j]=2
            
        for o = n_i+n_j>2 ? [cnt, randperm(n_i+n_j-2)+2] : cnt
                
            k=ixxs[o]
                
            # calculate the combined likelihood
                
            y_k = consts.datas[:,k]

            cstats,clikelihood = getlik(consts,priors,cstats,y_k,cns,true,true)
                
            clik=clik+clikelihood
            cns=cns+1

            cstats = update_Stats(cstats,y_k,cns,1)
                             
            # calculate individual set likelihoods
                
            if o.==1 || o.==2
#println(sstats[classhelp[o]])
                sstats[classhelp[o]],likelihood[o] = getlik(consts,priors,sstats[classhelp[o]],y_k,n_S[o]-1,true,true)
                setlik=setlik+likelihood[o]
                sstats[classhelp[o]] = update_Stats(sstats[classhelp[o]],y_k,n_S[o],1)
                continue
            end
            
            for ell = 1:2
                
                sstats[ell],likelihood[ell] = getlik(consts,priors,sstats[ell],y_k,n_S[ell],true,true)
                
            end
                
            likelihoods = exp(likelihood)

            m_S=copy(n_S)
            m_S[1]=c_i<=consts.sources ? m_S[1]-consts.ns[c_i] : m_S[1]
            m_S[2]=c_j<=consts.sources ? m_S[2]-consts.ns[c_j] : m_S[2]
            
            prob_i[1] = (m_S[1]*likelihoods[1])/sum(m_S'.*likelihoods) # the proba of choosing S_i for individual k
            prob_i[2] = 1-prob_i[1]
                
            setlik=setlik+likelihood[classhelp[o]]
            cprod=cprod*prob_i[classhelp[o]]
                
            n_S[classhelp[o]] = n_S[classhelp[o]]+1

            sstats[classhelp[o]] = update_Stats(sstats[classhelp[o]],y_k,n_S[classhelp[o]],1)
                                              
        end

        m_S=copy(n_S)
        m_S[1]=c_i<=consts.sources ? m_S[1]-consts.ns[c_i] : m_S[1]
        m_S[2]=c_j<=consts.sources ? m_S[2]-consts.ns[c_j] : m_S[2]
            
            
        M_H_prior = exp(lgamma(n_j+n_i)-(lgamma(m_S[1])+lgamma(m_S[2])))/alpha
        M_H_Lik =exp(clik-setlik)
        M_H_rat = M_H_prior*(M_H_Lik)*cprod
            
            if rand().<M_H_rat # accept ?
                #println("accept merge")
                # first update suff-stats of new merged group
                
                counts[merga] = n_i + n_j
                allcounts[merga] = allcounts[merga]+nmerge
               
                Stats[merga] = getlik(consts,priors,cstats,0,allcounts[merga],false,true)
        
                # then delete old table
                
                class_idd=classids_temp

                if antimerga>consts.sources
                    K_plus = K_plus-1
                
                    hits = [1:antimerga-1, antimerga+1:(K_plus+1)]
                
                counts[1:K_plus] = counts[hits]
                counts[K_plus+1] = 0

                allcounts[1:K_plus] = allcounts[hits]
                allcounts[K_plus+1] = 0
                
                Stats[1:K_plus] =  Stats[hits]
                Stats[K_plus+1] =  0

                else
                
                    Stats.means[:,antimerga] = consts.orgmeans[:,antimerga];
                    Stats.sum_squares[:,:,antimerga] = consts.orgsum_squares[:,:,antimerga];
                    counts[antimerga] = 0;
                    allcounts[antimerga] = consts.ns[antimerga];
                                       
                    Stats.log_det_cov[antimerga] = consts.orglog_det_cov[antimerga];
                    Stats.inv_cov[:,:,antimerga] = consts.orginv_cov[:,:,antimerga];
                    
                    K_plus = sum(allcounts.!=0);
                
                end
                
            end
            
    else # split - this is essentially the same thing only that the M-H is slightly different
            
            # reallocate, keeping track of the probas for the M-H ratio
            # for this I need the sufficient stats for the 'new' groups in
            # the merge case starting from one single obs in each group
            # [RANDOM SAMPLE THE ORDER! [Dahl 2003]]
        n_S=ones(Int16,2,1) # keeps track of the number of items at each 'new-old' component
           
        sstats = NORM(zeros(Float,D,2),zeros(Float,D,D,2),zeros(Float,D,D,2),zeros(Float,2))
        
        cstats = NORM( zeros(Float,D),zeros(Float,D,D),zeros(Float,D,D),zeros(Float))

        if c_i<=consts.sources
                    
            n_S[1] = n_S[1]+consts.ns[c_i]
            sstats.means[:,1] = consts.orgmeans[:,c_i] 
            sstats.sum_squares[:,:,1] = consts.orgsum_squares[:,:,c_i]

            cns = consts.ns[c_i]
            cstats.means = consts.orgmeans[:,c_i] 
            cstats.sum_squares = consts.orgsum_squares[:,:,c_i]
        else
            cns=0            
        end
       
           
            clik=0;cprod=1;setlik=0;likelihood=zeros(Float,(1,2))
            prob_i=zeros(Float,(2,1))
         
            classids_temp[ind[2]]=K_plus+1 # ind 1 marks new class
                 
            ixxs = find(class_idd.==c_j)
            ixxs=ixxs[ixxs.!=ind[1]]
            ixxs=ixxs[ixxs.!=ind[2]]
            ixxs=[ind[1], ind[2], ixxs]
            
            for o = n_j>2 ? [1, 2, randperm(n_j-2)+2] : [1, 2]
                
                k=ixxs[o]
                
                # calculate the combined likelihood
                y_k = consts.datas[:,k]
                
                cstats,clikelihood = getlik(consts,priors,cstats,y_k,cns,true,true)
                
                clik=clik+(clikelihood)
                cns=cns+1

                cstats = update_Stats(cstats,y_k,cns,1)
                
                # calculate individual set likelihoods

                m_S=deepcopy(n_S)
                if c_i<=consts.sources
                    m_S[1]=n_S[1]-consts.ns[c_i];
                end
                
                if o.==1 || o.==2                 
                    
                    sstats[o],likelihood[o] = getlik(consts,priors,sstats[o],y_k,n_S[o]-1,true,true)
                    sstats[o] = update_Stats(sstats[o],y_k,n_S[o],1)
                    setlik=setlik+likelihood[o]
                    continue                   
                end
                
                for ell = 1:2

                    sstats[ell],likelihood[ell] = getlik(consts,priors,sstats[ell],y_k,n_S[ell],true,true)
                                    
                end
                
                likelihoods = exp(likelihood)
           
                
                prob_i[1] = (m_S[1]*likelihoods[1])/sum(m_S'.*likelihoods) # the proba of choosing S_i for individual k
                prob_i[2] = 1-prob_i[1]
                
                if rand().<prob_i[1] # S_i is chosen, this time I need to keep track of which one k was allocated to                    
                    setlik=setlik+likelihood[1]
                    cprod=cprod*prob_i[1]
                    
                    
                    n_S[1] = n_S[1]+1
                    
                    sstats[1] = update_Stats(sstats[1],y_k,n_S[1],1)
                    
                else # S_j is chosen
                    
                    setlik=setlik+likelihood[2]
                    cprod=cprod*prob_i[2]
                    classids_temp[k]=K_plus+1
                    n_S[2] = n_S[2]+1
                    
                    sstats[2] = update_Stats(sstats[2],y_k,n_S[2],1)
                    
                end
                
            end

        m_S=deepcopy(n_S)
        if c_i<=consts.sources
            m_S[1]=n_S[1]-consts.ns[c_i];
        end
            
        M_H_prior = exp((lgamma(m_S[1])+lgamma(m_S[2]))-lgamma(n_j))*alpha
        M_H_Lik = exp(setlik-clik)
        M_H_rat = M_H_prior*(M_H_Lik)*(1/cprod)
            
        if rand().<M_H_rat #&& any(classids_temp!=class_idd) # accept ?
                # println("accept split")
                # first update suff-stats of new groups
                if c_i<=consts.sources
                    counts[c_j] =n_S[1]-consts.ns[c_i];
                else                
                    counts[c_j] = n_S[1]
                end
                
                counts[K_plus+1] = n_S[2]
                allcounts[c_j] = n_S[1]
                allcounts[K_plus+1] = n_S[2]

                Stats[c_j] = getlik(consts,priors,sstats[1],0,allcounts[c_j],false,true)
                Stats[K_plus+1] = getlik(consts,priors,sstats[2],0,allcounts[K_plus+1],false,true)
                
                K_plus=K_plus+1
                class_idd=classids_temp
            end
            
        end

       return(class_idd,K_plus,Stats,counts,allcounts)
        
end
