using Distributions

function preclass(consts,Stats,priors,allcounts,counts,class_id)

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
               
        Stats.means[:,choose] = Stats.means[:,choose] + (1/n)*(consts.datas[:,i]-Stats.means[:,choose]);
        Stats.sum_squares[:,:,choose] = Stats.sum_squares[:,:,choose] + consts.datas[:,i]*consts.datas[:,i]'
        
        
        k_n = priors.k_0+n;
        v_n = priors.v_0+n;
        
        zm_Y = Stats.means[:,choose]-priors.mu_0
        SS = Stats.sum_squares[:,:,choose]-n*(Stats.means[:,choose]*Stats.means[:,choose]')
        lambda_n = priors.lambda_0 + SS + priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    
        Stats.log_det_cov[choose] = log(det(Sigma))
        Stats.inv_cov[:,:,choose] = inv(Sigma)
    
    end

    return (Stats,allcounts,counts)
end

# pre-allocate stats and consts


function preallocate(consts,Stats,priors,allcounts)

    is=unique(consts.labels,true);

   for i=1:consts.sources
    
        Stats.means[:,i] = mean(consts.baseline[:,consts.labels.==is[i]],2)
        consts.orgmeans[:,i] =  Stats.means[:,i]
        allcounts[i] = size(consts.baseline[:,consts.labels.==is[i]],2)
        consts.ns[i] = size(consts.baseline[:,consts.labels.==is[i]],2)
        consts.orgsum_squares[:,:,i] = (baseline[:,consts.labels.==is[i]]-repmat(mean(consts.baseline[:,consts.labels.==is[i]],2),1,allcounts[i]))*(consts.baseline[:,consts.labels.==is[i]]-repmat(mean(consts.baseline[:,consts.labels.==is[i]],2),1,allcounts[i]))'
        Stats.sum_squares[:,:,i] =  consts.orgsum_squares[:,:,i]
    
        consts.orgsum_squares[:,:,i] = consts.orgsum_squares[:,:,i]+consts.ns[i]*(Stats.means[:,i]*Stats.means[:,i]');
        Stats.sum_squares[:,:,i] = Stats.sum_squares[:,:,i]+consts.ns[i]*(Stats.means[:,i]*Stats.means[:,i]');
    
        n = size(consts.baseline[:,consts.labels.==is[i]],2)
        k_n = priors.k_0+n
        v_n = priors.v_0+n
    
        zm_Y = Stats.means[:,i]-priors.mu_0
        SS = Stats.sum_squares[:,:,i]-n*(Stats.means[:,i]*Stats.means[:,i]')
        lambda_n = priors.lambda_0 + SS +  priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
        Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    
        Stats.log_det_cov[i] = log(det(Sigma))
        consts.orglog_det_cov[i] =  Stats.log_det_cov[i]
        Stats.inv_cov[:,:,i] = inv(Sigma)
        consts.orginv_cov[:,:,i] =  Stats.inv_cov[:,:,i]
    end

    return(consts,Stats,allcounts)

end

#unique

function unique{T}(A::AbstractArray{T}, sorted::Bool)
    dd = Dict{T, Bool}()
    for a in A dd[a] = true end
    sorted? sort!(keys(dd)): keys(dd)
end

# update statistics Normal case

function update_Stats(Stats::NORM,y,counts,m)

    if m<0
        Stats.means = (1/counts)*((counts+1)*Stats.means- y);
        Stats.sum_squares = Stats.sum_squares - y*y'
    else
        Stats.means=  Stats.means+ (1/counts)*(y- Stats.means);
        Stats.sum_squares =  Stats.sum_squares + y*y'
    end
    return Stats
end


# student log lik
function getlik(consts::STUD,priors::MNIW,Stats::NORM,y,n,lik::Bool,suffs::Bool)
         
    m_Y = Stats.means
    mu = priors.k_0/(priors.k_0+n)*priors.mu_0 + n/(priors.k_0+n)*m_Y
    k_n = priors.k_0+n
    v_n = priors.v_0+n

    S = (Stats.sum_squares - n*m_Y*m_Y')
    zm_Y = m_Y-priors.mu_0
    lambda_n = priors.lambda_0 + S  +  priors.k_0*n/(priors.k_0+n)*zm_Y*zm_Y'

    Sigma = lambda_n*(k_n+1)/(k_n*(v_n-consts.D+1))
    v = v_n-consts.D+1
   
    vd = v+consts.D
    d2 = consts.D/2

    if suffs

        Stats.log_det_cov =  log(det(Sigma))
        Stats.inv_cov  = inv(Sigma)

        if lik
        
            lp = (consts.pc_gammaln_by_2[vd] - (consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] + d2*consts.pc_log_pi) - .5*Stats.log_det_cov-(vd/2)*log(1+(1/v)*(y-mu)'*Stats.inv_cov*(y-mu)))[1]

            return Stats,lp
        else
            return Stats
        end
    else
        lp = (consts.pc_gammaln_by_2[vd] - (consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] + d2*consts.pc_log_pi) - .5*Stats.log_det_cov-(vd/2)*log(1+(1/v)*(y-mu)'*Stats.inv_cov*(y-mu)))[1]
    
        return lp
    end
    
end


### p under prior alone

function p_for_1(consts::STUD,priors::MNIW,N,datas,p_under_prior_alone)

 Sigma = (priors.lambda_0*(priors.k_0+1)/(priors.k_0*(priors.v_0-consts.D+1)))'
    v = priors.v_0-consts.D+1
    mu = priors.mu_0
    log_det_Sigma = log(det(Sigma))
    inv_Sigma = Sigma^-1
    vd = v+consts.D
    d2=consts.D/2
    for i=1:N
        y = datas[:,i]
              
        lp =consts.pc_gammaln_by_2[vd] - ( consts.pc_gammaln_by_2[v] + d2*consts.pc_log[v] +  d2*consts.pc_log_pi) -.5*log_det_Sigma- (vd/2)*log(1+(1/v)*(y-mu)'*inv_Sigma*(y-mu))
        
        p_under_prior_alone[i] = lp[1]

    end

    return p_under_prior_alone
end

## display time remaining

function disptime(total_time,time_1_iter,iter,thin,num_iters,K_record)

total_time = total_time + time_1_iter
    if iter.==1
       println(string("Iter: ",dec(iter),"/",dec(num_iters)))
    elseif mod(iter,thin*5).==0
        E_K_plus = mean(K_record[1:int(iter/thin)])
        rem_time = (time_1_iter*.05 + 0.95*(total_time/iter))*num_iters-total_time
        if rem_time < 0
            rem_time = 0
        end
        println(string("Iter: ",dec(iter),'/',dec(num_iters),", Rem. Time: ", secs2hmsstr(rem_time),", mean[K^+]",string(E_K_plus)))
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

 g_alpha = rand(Gamma(alpha+1))/2
 g_N = rand(Gamma(N))/2
 nu = g_alpha/(g_alpha+g_N)
 pis=(a_0+K_plus-1)/(a_0+K_plus-1+N*(b_0-log(nu)))
 alpha = (pis*(rand(Gamma(a_0+K_plus))/(b_0-log(nu))))+((1-pis)*(rand(Gamma(a_0+K_plus-1))/(b_0-log(nu))))

 end

 # update prior

 function update_prior(consts::STUD,K_plus,stats::NORM,counts,priors::MNIW)

     muu = zeros(Float,(D,K_plus))
     sums=0
     invsig = zeros(Float,(D,D,K_plus))
     sumsig=0
     for k =1:K_plus
                
         m_Y = stats.means[:,k]
         n=counts[k]
         mu_n = priors.k_0/(priors.k_0+n)*priors.mu_0 + n/(priors.k_0+n)*m_Y
         SS = (stats.sum_squares[:,:,k] - n*(m_Y*m_Y'))
         zm_Y = m_Y-priors.mu_0
         lambda_n = priors.lambda_0 + SS + priors.k_0*n/(priors.k_0+n)*(zm_Y)*(zm_Y)'
                
         v_n=priors.v_0+n

         #simualte from wishard
         for i = 1:v_n
             try
                 A = chol(lambda_n)*randn(consts.D)
             catch
                 return(priors.k_0,priors.mu_0)
             end
         
             
             invsig[:,:,k] += A*A'
         end
                
         # simulate from mvnorm
         A = chol(inv(invsig[:,:,k])/(priors.k_0+n))
         zs = randn(consts.D)
         
         muu[:,k]=mu_n+A*zs              
         #println( muu[:,k])
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
         
         sums += (muu[:,k]-mu_0)'*invsig[:,:,k]*(muu[:,k]-mu_0)
     end
     
     k_0 = rand(Gamma((K_plus+1)/2))*((1/(sums)+1)/2)[1]
    
     return(k_0,mu_0[:,1])

 end
