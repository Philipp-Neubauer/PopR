# julia code main wrapper

bl=int(ARGS[1])
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


    # define normal set types

     
    if bl.==1
    
        @everywhere require("fixtype.jl")

    else


        @everywhere require("define_types.jl")
      
    end
 

#################################################
######### --- set up outputs ---- ###############
#################################################
    
np = nprocs()
nit=int((np-1)*numiters/thin);

class_ids=Array(Float,(size(datas,2),nit))
k_0s=Array(Float,nit)
K_record=Array(Int8,nit)
alpha_record=Array(Float,nit)

#################################################
######### --- RUN IT ----########################
#################################################

if bl.==1

    @everywhere require("DPM_sampler_fix.jl")
    
    inp = {numiters,thin,Stats,priors,consts}
    M = {inp for i=1:np-1}
        
    outs = pmap(DPM_sampler_fix,M)
   
else

    @everywhere require("DPM_sampler.jl")
    
    inp = {datas,numiters,thin,Stats,priors,consts};
    M = {inp for i=1:np-1};
        
    outs = pmap(DPM_sampler,M)
   
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
