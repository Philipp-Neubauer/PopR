function  MM_sampler(num_iters,thin,cond,stats,priors,consts)

    
    
    require(string(CP,"/Julia_code_support.jl"))

    # initialize structures

    nit=int(num_iters/thin)
   
    counts = zeros(Int,consts.sources)
    allcounts = zeros(Int,consts.sources,1)
    class_id = Array(Int,consts.N)
    prop = ones(Float,consts.sources)
    class_ids = Array(Int,(consts.N,nit))
    props = Array(Float,(consts.sources,nit))

    predlik=ones(Float,(consts.N,consts.sources))'
    proba=Array(Float,(consts.N,consts.sources))
    probas=Array(Float,(consts.N,consts.sources,nit))
    
    consts,stats,allcounts = preallocate(consts,stats,priors,allcounts)
    stats,allcounts,counts = preclass(consts,stats,priors,allcounts,counts,class_id)
   
    tic()
    totaltime=0

 # get likelihood for each ind if cond

    #if cond ==1
 
        for i=1:consts.N
            for j=1:consts.sources

                predlik[j,i] = getlik(consts,priors,stats[j],consts.datas[:,i],consts.ns[j],true,false)

            end
        end

    #end
    
    ## start MCMC ---- 
    for iter=1:num_iters
    #println(iter)
        for i=1:consts.N

           counts[class_id[i]] = counts[class_id[i]] - 1
        
            
            for j=1:consts.sources      
                # get likelihood for each ind if uncond
                if cond == 0

                    # take out in from stats if uncond
                    if class_id[i] == j
                        allcounts[j] = allcounts[j]-1
                        stats[j] = update_stats(stats[j],consts.datas[:,i],allcounts[j],-1)
                    end
                    
                    (dum,predlik[j,i]) = getlik(consts,priors,stats[j],consts.datas[:,i],allcounts[j],true,true)
                  
                end
                               
            end
            
             # allocate
                
            proba[i,:] = exp(predlik[:,i]).*prop/sum(exp(predlik[:,i]).*prop)
            prob = cumsum(proba[i,:]')
            class_id[i] = find(prob.>rand())[1]
            counts[class_id[i]] = counts[class_id[i]] + 1    
            
            # add to stats if uncond
            if cond == 0
                allcounts[class_id[i]] = allcounts[class_id[i]] + 1    
                stats[class_id[i]]=update_stats(stats[class_id[i]],consts.datas[:,i],allcounts[class_id[i]],1)
            end
        end 
       

       # update source proportions/priors

       prop = rand(Dirichlet((counts+1/consts.sources)))
       
        # save parameter values
        
        if mod(iter,thin)==0
           class_ids[:,int(iter/thin)]=class_id
           props[:,int(iter/thin)]=prop
           probas[:,:,int(iter/thin)]=proba
        end
  
        # timer

        time_1_iter = toq();
        #println(prop)
        totaltime = disptime(totaltime,time_1_iter,iter,thin,num_iters)
        tic()
        
    end

    return(class_ids,props,reshape(mean(probas,3),consts.N,consts.sources))

end
