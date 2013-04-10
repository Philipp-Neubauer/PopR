# julia code main wrapper

bl=int(ARGS[1])
thin =int(ARGS[3])
numiters =int(ARGS[2])

Typeof=ARGS[4]

cd(ARGS[5])


#d( "/home/philbert/Papersampler")

# define type alias for floats same as for Int

#if is(Int,Int64)
    typealias Float Float64
#else
#    typealias Float Float32
#end


datas =readdlm("datas.csv",",",Float)
global datas=datas[2:end,2:end]'

(D, N) = size(datas)

########################################
# define structures for normal model####
########################################


#if(Typeof=="N")

    single_priors = readdlm("single_priors.csv",",",Float)
    single_priors=single_priors[2:end,2]
    
    matrix_priors = readdlm("matrix_priors.csv",",",Float)
    matrix_priors=matrix_priors[2:end,2:end]
   
    # define normal set types

    pc_max_ind=int(1e5)
    
    if bl.==1
    
        require("fixtype.jl")

    else


require("define_types.jl")
      
    end
 

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
  
   @everywhere require("DPM_sampler_fix.jl")
   #@everywhere load("define_types.jl")
   for n=1:np

       push!(outs,fetch(@spawn DPM_sampler_fix(numiters,thin,Stats,priors,consts)))
 
   end
   
else

   @everywhere require("DPM_sampler.jl")
   # @everywhere require("define_types.jl")
   for n=1:np
    
       push!(outs,fetch(@spawn DPM_sampler(datas,numiters,thin,Stats,priors,consts)))

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

writecsv("source_ids.csv",class_ids)
writecsv("K_record.csv",K_record)
writecsv("gammas.csv",alpha_record)
writecsv("k_0s.csv",k_0s)
