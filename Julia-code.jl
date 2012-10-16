function call_DPM(bl,num_iters,thin,wd)

# julia code main wrapper function

cd($wd)

single_priors = csvread("priors_singles.csv")
matrix_priors = csvread("priors_matrix.csv")
datas = csvread("datas.csv")

outs={}
np = nprocs()

class_ids=Array(Float32,np*numiters/thin,size(datas,1))
k_0s=Array(Float32,np*numiters/thin)
K_record=Array(Int8,np*numiters/thin)
alpha_record=Array(Float32,np*numiters/thin)

if bl.==1
   baseline = csvread("baseline.csv")
   @everywhere load("DPM_sampler_fix.jl")
   for n=1:np
       push(outs,@spawn DPM_sampler_fix(datas,baseline_single_priors,matrix_priors,num_iters,thin)
 
   end
else
   @everywhere load("DPM_sampler.jl")
   for n=1:np
    
       push(outs,@spawn DPM_sampler(datas,single_priors,matrix_priors,num_iters,thin)

   end
end

 while length(outs)>0
     res = fetch(pop(outs))
     class_ids[((n-1)*numiters/thin+1):((n-1)*numiters/thin+numiters/thin),:] = res[1]
     k_0s[((n-1)*numiters/thin+1):((n-1)*numiters/thin+numiters/thin)] = res[2]
     K_record[((n-1)*numiters/thin+1):((n-1)*numiters/thin+numiters/thin)] =res[3]
     alpha_record[((n-1)*numiters/thin+1):((n-1)*numiters/thin+numiters/thin)] = res[4]
end

# write out simualtions

csv_write("source_ids.csv",class_ids)
csv_write("K_record.csv",K_record)
csv_write("gammas.csv",alpha_record)
csv_write("k_0s.csv",k_0s)

end

outs={}
for n=1:2
   push(outs, @spawn rand2(2,2))
   end
   
 while length(outs)>0
 res = fetch(pop(outs))
 end
   
