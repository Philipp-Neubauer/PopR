# julia code main wrapper
using Distributions
condi=int(ARGS[1])
thin =int(ARGS[3])
numiters =int(ARGS[2])

Typeof=ARGS[4]

cd(ARGS[5])


@everywhere typealias Float Float64

@everywhere datas=readdlm("datas.csv",",",Float)[2:end,2:end]'

@everywhere single_priors = readdlm("single_priors.csv",",",Float)[2:end,2]
    
@everywhere matrix_priors = readdlm("matrix_priors.csv",",",Float)[2:end,2:end]


@everywhere const pc_max_ind=int(1e5)
@everywhere const (D, N) = size(datas)
    
@everywhere    require("fixtype.jl")   
 

#################################################
######### --- set up outputs ---- ###############
#################################################
   nit=int(numiters/thin);

class_ids=Array(Float,(size(datas,2),nit))
props=Array(Float,(consts.sources,nit))
probas=Array(Float,(size(datas,2),consts.sources))


#################################################
######### --- RUN IT ----########################
#################################################

@everywhere require("MM_sampler.jl")
#@everywhere load("define_types.jl")


    outs = MM_sampler(numiters,thin,condi,Stats,priors,consts)
 

# get outputs

nit=int(numiters/thin);
n=1
  
        class_ids[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[1]
        props[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[2]
        probas = outs[3]
  

#################################################
######### --- Write out ----#####################
#################################################

writecsv("source_ids.csv",class_ids)
writecsv("proportions.csv",props)
writecsv("post_probas.csv",probas)
