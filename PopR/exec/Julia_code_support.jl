using Distributions

function preclass(consts,stats,priors,allcounts,counts,class_id)

    dists=Array(Float,consts.sources)
    for i=1:consts.N
        
        for j=1:(consts.sources)
            dists[j] =  sqrt(sum((consts.datas[:,i]-consts.orgmeans[:,j]).^2));
        end
    
        choose=find(rand().>=(1-cumsum(dists./sum(dists))))[1]
        allcounts[choose]=allcounts[choose]+1
        counts[choose]=counts[choose]+1
        class_id[i]=choose;
        
        n = allcounts[choose]
               
        stats.means[:,choose] = stats.means[:,choose] + (1/n)*(consts.datas[:,i]-stats.means[:,choose]);
        stats.sum_squares[:,:,choose] = stats.sum_squares[:,:,choose] + consts.datas[:,i]*consts.datas[:,i]'
        
        
        k_n = priors.k_0+n;
        v_n = priors.v_0+n;
        
        zm_Y = stats.means[:,choose]-priors.mu_0
        SS = stats.sum_squares[:,:,choose]-n*(stats.means[:,choose]*stats.means[:,choose]')
        lambda_n = priors.lambda_0 + SS + priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    
        stats.inv_cov[:,:,choose] = cholfact(Sigma)[:U]
      
        stats.log_det_cov[choose] = sum(log(diag(stats[choose].inv_cov)))    
    end

    return (stats,allcounts,counts)
end

# pre-allocate stats and consts


function preallocate(consts,stats,priors,allcounts)

    is=sort(unique(consts.labels));

   for i=1:consts.sources
    
        stats.means[:,i] = mean(consts.baseline[:,consts.labels.==is[i]],2)
        consts.orgmeans[:,i] =  stats.means[:,i]
        allcounts[i] = size(consts.baseline[:,consts.labels.==is[i]],2)
        consts.ns[i] = size(consts.baseline[:,consts.labels.==is[i]],2)
        consts.orgsum_squares[:,:,i] = (baseline[:,consts.labels.==is[i]]-repmat(mean(consts.baseline[:,consts.labels.==is[i]],2),1,allcounts[i]))*(consts.baseline[:,consts.labels.==is[i]]-repmat(mean(consts.baseline[:,consts.labels.==is[i]],2),1,allcounts[i]))'
        stats.sum_squares[:,:,i] =  consts.orgsum_squares[:,:,i]
    
        consts.orgsum_squares[:,:,i] = consts.orgsum_squares[:,:,i]+consts.ns[i]*(stats.means[:,i]*stats.means[:,i]');
        stats.sum_squares[:,:,i] = stats.sum_squares[:,:,i]+consts.ns[i]*(stats.means[:,i]*stats.means[:,i]');
    
        n = size(consts.baseline[:,consts.labels.==is[i]],2)
        k_n = priors.k_0+n
        v_n = priors.v_0+n
   
        zm_Y = stats.means[:,i]-priors.mu_0
        SS = stats.sum_squares[:,:,i]-n*(stats.means[:,i]*stats.means[:,i]')
        lambda_n = priors.lambda_0 + SS +  priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    
        stats.inv_cov[:,:,i] = cholfact(Sigma)[:U]
        consts.orginv_cov[:,:,i] =  stats.inv_cov[:,:,i]
        stats.log_det_cov[i] = sum(log(diag(stats[i].inv_cov)))
        consts.orglog_det_cov[i] =  stats.log_det_cov[i]
    end

    return(consts,stats,allcounts)

end



# update statistics Normal case

function update_stats(stats::NORM,y,counts,m)

    if m<0
        stats.means = (1/counts)*((counts+1)*stats.means- y);
        stats.sum_squares = stats.sum_squares - y*y'
    else
        stats.means=  stats.means+ (1/counts)*(y- stats.means);
        stats.sum_squares =  stats.sum_squares + y*y'
    end
    return stats
end


# student log lik
function getlik(consts::STUD,priors::MNIW,stats::NORM,y,n,lik::Bool,suffs::Bool)
         
    m_Y = stats.means
    mu = priors.k_0/(priors.k_0+n)*priors.mu_0 + n/(priors.k_0+n)*stats.means
    k_n = priors.k_0+n
    v_n = priors.v_0+n

    S = (stats.sum_squares - n* stats.means*stats.means')
    zm_Y = m_Y-priors.mu_0
    lambda_n = priors.lambda_0 + S  +  priors.k_0*n/(priors.k_0+n)*zm_Y*zm_Y'

    Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    v = v_n-consts.D+1
   
    vd = v+consts.D
    d2 = consts.D/2

#println(Sigma)
    
    if suffs

        stats.inv_cov  = cholfact(Sigma)[:U]
        stats.log_det_cov = sum(log(diag(stats.inv_cov)))
        
        if lik
            u = y-mu
#println(u)
            
            z = Cholesky(stats.inv_cov,'U') \ u  # This is equivalent to inv(cov) * u, but much faster
            #println(z)
           # print(z,'\n')
            lp = consts.pc_gammaln_by_2[vd] - (consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] + d2*consts.pc_log_pi) - stats.log_det_cov-(vd/2)*log(1+(1/v)*dot(u,z))
                  
            return stats,lp
        else
            return stats
        end
    else

         u = y-mu
         z = Cholesky(stats.inv_cov,'U') \ u  # This is equivalent to inv(cov) * u, but much faster
                 
        lp = consts.pc_gammaln_by_2[vd] - (consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] + d2*consts.pc_log_pi) - stats.log_det_cov-(vd/2)*log(1+(1/v)*dot(u,z))
    
        return lp
    end
    
end


### p under prior alone

function p_for_1(consts::STUD,priors::MNIW,N,datas,p_under_prior_alone)

 Sigma = (priors.lambda_0*(priors.k_0+1)/(priors.k_0*(priors.v_0-consts.D+1)))'
    v = priors.v_0-consts.D+1
    mu = priors.mu_0
    #println(Sigma)
 inv_Sigma = cholfact(Sigma)
    log_det_Sigma =sum(log(diag( inv_Sigma[:U])))
   
    vd = v+consts.D
    d2=consts.D/2
    for i=1:N
        y = datas[:,i]
               u = y-mu
z = inv_Sigma \ u
          lp = consts.pc_gammaln_by_2[vd] - (consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] + d2*consts.pc_log_pi) - log_det_Sigma-(vd/2)*log(1+(1/v)*dot(u,z))
        
        p_under_prior_alone[i] = lp

    end

    return p_under_prior_alone
end

## display time remaining

function disptime(total_time,time_1_iter,iter,thin,num_iters,K_record)

total_time = total_time + time_1_iter
    if iter.==1
       println(string("Iter: ",dec(iter),"/",dec(num_iters)))
    elseif mod(iter,thin*100).==0
        E_K_plus = round(mean(K_record[1:int(iter/thin)]),2)
        rem_time = (time_1_iter*.05 + 0.95*(total_time/iter))*num_iters-total_time
        if rem_time < 0
            rem_time = 0
        end
        println(string("Iter: ",dec(iter),'/',dec(num_iters),", Rem. Time: ", secs2hmsstr(rem_time),", mean[K^+]",string(E_K_plus)))
    end

    return(total_time)

end

## display time remaining

function disptime(total_time,time_1_iter,iter,thin,num_iters)

total_time = total_time + time_1_iter
    if iter.==1
       println(string("Iter: ",dec(iter),"/",dec(num_iters)))
    elseif mod(iter,thin*100).==0
        rem_time = (time_1_iter*.05 + 0.95*(total_time/iter))*num_iters-total_time
        if rem_time < 0
            rem_time = 0
        end
        println(string("Iter: ",dec(iter),'/',dec(num_iters),", Rem. Time: ", secs2hmsstr(rem_time)))
    end

    return(total_time)

end

# convert secs into dhminsec

function secs2hmsstr(secs)

days = ifloor(secs/(3600*24))
rem = (secs-(days*3600*24))
hours = ifloor(rem/3600)
rem = rem -(hours*3600)
minutes = ifloor(rem/60)
rem = rem - minutes*60
secs = ifloor(rem)

    if days .== 0
        str =  string(dec(hours),"h:",dec(minutes),"min:",dec(secs),"sec")
   elseif days .==1
        str = string( "1 Day + ",dec(hours),"h:",dec(minutes),"min:",dec(secs),"sec")
    else
        str = string(dec(days),"days:",dec(hours),"h:",dec(minutes),"min:",dec(secs),"sec")
end

end


# update alpha

function update_alpha(alpha,N,K_plus,a_0,b_0)

 nu = rand(Beta(alpha+1,N))
 pis=(a_0+K_plus-1)/(a_0+K_plus-1+N*(b_0-log(nu)))
 alpha = (pis*(rand(Gamma(a_0+K_plus))/(b_0-log(nu))))+((1-pis)*(rand(Gamma(a_0+K_plus-1))/(b_0-log(nu))))

 end

 # update prior

 function update_prior(consts::STUD,K_plus,stats::NORM,counts,priors::MNIW)

     muu = zeros(Float,(D,K_plus))
     sums=0
     sig = zeros(Float,(D,D,K_plus))
     invsig= zeros(Float,(D,D,K_plus))
     sumsig=0
     #print(priors.k_0,'\n')
     for k =1:K_plus
                
       
         n=counts[k]
         mu_n = priors.k_0/(priors.k_0+n)*priors.mu_0 + n/(priors.k_0+n)*stats.means[:,k]
         SS = (stats.sum_squares[:,:,k] - n*(stats.means[:,k]*stats.means[:,k]'))
         zm_Y = stats.means[:,k]-priors.mu_0
         lambda_n = priors.lambda_0 + SS + priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
                
         v_n=priors.v_0+n

         try
             sig[:,:,k] = rand(InverseWishart(v_n,lambda_n))

         catch
             return(priors.k_0,priors.mu_0)
         end
     
         # simulate from mvnorm
         A = chol(sig[:,:,k]/(priors.k_0+n))
         zs = randn(consts.D)
         
         muu[:,k]=mu_n+A*zs              
         #println( muu[:,k])
         invsig[:,:,k] = inv(sig[:,:,k])
         sums += invsig[:,:,k]*muu[:,k]
         sumsig += invsig[:,:,k]
     end
            
     meansig=inv(sumsig)
     
     A = chol(meansig./priors.k_0)
     zs = randn(consts.D)
     
     mu_0=(meansig/sums')+A*zs    
     # println( mu_0)
     
     sums=0
     for k=1:K_plus
         
         sums += (muu[:,k]-mu_0)'*(invsig[:,:,k])*(muu[:,k]-mu_0)
     end

     #print(inv(invsig[:,:,1]),'\n')
     
     k_0 = rand(Gamma((K_plus+priors.ak_0)/2))*((sums+priors.bk_0)/2)[1]
    #println(k_0)
     return(k_0,mu_0[:,1])

 end


