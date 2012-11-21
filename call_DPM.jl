# julia code main wrapper

bl=int(ARGS[1])
thin =int(ARGS[3])
numiters =int(ARGS[2])

Typeof=ARGS[4]

cd(ARGS[5])


#cd( "C:/Users/Philipp/Work/Otolith_Code/HBM/DPMpaper/Papersampler")

# define type alias for floats same as for Int

if is(Int,Int64)
    typealias Float Float64
else
    typealias Float Float32
end

datas =dlmread("datas.csv",",",Float)
global datas=datas[2:end,2:end]'

(D, N) = size(datas)

########################################
# define structures for normal model####
########################################


#if(Typeof=="N")

    single_priors = dlmread("single_priors.csv",",",Float)
    single_priors=single_priors[2:end,2]
    
    matrix_priors = dlmread("matrix_priors.csv",",",Float)
    matrix_priors=matrix_priors[2:end,2:end]
   
    # define normal set types

    pc_max_ind=int(1e5)
    
    if bl.==1
    
        load("fixtype.jl")

    else

        load("define_types.jl")
      
    end
  
     
    # prior composite type
    type MNIW
        
        a_0::Float
        b_0::Float
        k_0::Float
        v_0::Float
        mu_0::Array{Float,1}
        lambda_0::Array{Float,2}
        
    end
    
    priors= MNIW(single_priors[1],single_priors[2],single_priors[3],single_priors[4],single_priors[5:end],matrix_priors)
    
    # stats composite type
    type NORM
        
        means
        sum_squares
        inv_cov
        log_det_cov
        
    end
    
     # set up NORM type stats
    Stats=NORM(zeros(Float,(D,N)),zeros(Float,(D,D,N)),Array(Float,(D,D,N)),Array(Float,(N)))
    
    # define how to access subsets of individuals
    
    function ref(A::NORM,k::Any)

        NORM(A.means[:,k],A.sum_squares[:,:,k],A.inv_cov[:,:,k],A.log_det_cov[k])
            
    end

     # define how to assign subsets of individuals
    
    function assign(A::NORM,B::NORM,k::Any)

        A.means[:,k]=B.means
        A.sum_squares[:,:,k] = B.sum_squares
        A.inv_cov[:,:,k] = B.inv_cov
        A.log_det_cov[k] = B.log_det_cov
    
    end

    function assign(A::NORM,B::Number,k::Any)

        A.means[:,k]=B
        A.sum_squares[:,:,k] = B
        A.inv_cov[:,:,k] = B
        A.log_det_cov[k] = B
    
    end
    
#end

#################################################
######### --- set up outputs ---- ###############
#################################################
    
outs={}
np = nprocs()
nit=int(np*numiters/thin);

class_ids=Array(Float,(size(datas,2),nit))
k_0s=Array(Float,nit)
K_record=Array(Int8,nit)
alpha_record=Array(Float,nit)

#################################################
######### --- RUN IT ----########################
#################################################

if bl.==1
  
   @everywhere load("DPM_sampler_fix.jl")
   #@everywhere load("define_types.jl")
   for n=1:np

       push(outs,fetch(@spawn DPM_sampler_fix(numiters,thin,Stats,priors,consts)))
 
   end
   
else

   @everywhere load("DPM_sampler.jl")
   # @everywhere load("define_types.jl")
   for n=1:np
    
       push(outs,fetch(@spawn DPM_sampler(datas,numiters,thin,Stats,priors,consts)))

   end
   
end

# get outputs

nit=int(numiters/thin);
n=0
    for i=1:length(outs)
        n=n+1
        class_ids[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][1]
        k_0s[((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][2]
        K_record[((n-1)*nit+1):((n-1)*nit+nit)] =outs[i][3]
        alpha_record[((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][4]
    end

#################################################
######### --- Write out ----#####################
#################################################

csvwrite("source_ids.csv",class_ids)
csvwrite("K_record.csv",K_record)
csvwrite("gammas.csv",alpha_record)
csvwrite("k_0s.csv",k_0s)
