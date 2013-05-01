# julia code main wrapper
using Distributions
cond=int(ARGS[1])
thin =int(ARGS[3])
numiters =int(ARGS[2])

Typeof=ARGS[4]

cd(ARGS[5])


#cd( "/home/philbert/Papersampler")

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
    
    require("fixtype.jl")   
 

#################################################
######### --- set up outputs ---- ###############
#################################################
    
outs={}
np = nprocs()
nit=int(np*numiters/thin);

class_ids=Array(Float,(size(datas,2),nit))
props=Array(Float,(consts.sources,nit))
probas=Array(Float,(size(datas,2),consts.sources))


#################################################
######### --- RUN IT ----########################
#################################################

@everywhere require("MM_sampler.jl")
#@everywhere load("define_types.jl")
for n=1:np

    push!(outs,fetch(@spawn MM_sampler(numiters,thin,cond,Stats,priors,consts)))
 
end


# get outputs

nit=int(numiters/thin);
n=0
    for i=1:length(outs)
        n=n+1
        class_ids[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][1]
        props[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][2]
        probas = outs[i][3]
    end

#################################################
######### --- Write out ----#####################
#################################################

writecsv("source_ids.csv",class_ids)
writecsv("proportions.csv",props)
writecsv("post_probas.csv",probas)
