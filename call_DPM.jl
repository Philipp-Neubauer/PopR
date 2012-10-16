# julia code main wrapper

bl=int(ARGS[1])
thin = int(ARGS[3])
numiters = int(ARGS[2])

cd(ARGS[4])

load("DPM_sampler.jl")


datas =dlmread("datas.csv",",",Float64)
datas=datas[2:end,2:end]'

outs={}
np = nprocs()
nit=int(np*numiters/thin);

class_ids=Array(Float64,(size(datas,2),nit))
k_0s=Array(Float64,nit)
K_record=Array(Int8,nit)
alpha_record=Array(Float64,nit)


if bl.==1
  
   @everywhere load("DPM_sampler_fix.jl")
   for n=1:np
       push(outs,fetch(@spawn DPM_sampler_fix(numiters,thin)))
 
   end
else
   @everywhere load("DPM_sampler.jl")
   for n=1:np
    
       push(outs,fetch(@spawn DPM_sampler(numiters,thin)))

   end
end

nit=int(numiters/thin);
n=0
    for i=1:length(outs)
        n=n+1
        class_ids[:,((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][1]
        k_0s[((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][2]
        K_record[((n-1)*nit+1):((n-1)*nit+nit)] =outs[i][3]
        alpha_record[((n-1)*nit+1):((n-1)*nit+nit)] = outs[i][4]
    end

# write out simualtions

csvwrite("source_ids.csv",class_ids)
csvwrite("K_record.csv",K_record)
csvwrite("gammas.csv",alpha_record)
csvwrite("k_0s.csv",k_0s)
