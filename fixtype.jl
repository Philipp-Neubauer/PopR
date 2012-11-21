  baseline =dlmread("baseline.csv",",",Float)
        global baseline=baseline[2:end,2:end]'

        label = dlmread("labels.csv",",",Float)
        label = label[2:end,2]
        label=int(label)  
        

function unique{T}(A::AbstractArray{T}, sorted::Bool)
    dd = Dict{T, Bool}()
    for a in A dd[a] = true end
    sorted? sort!(keys(dd)): keys(dd)
end

        uniqs=unique(label,true)
        labels=Array(Int,length(label))
        for i=1:length(uniqs)
            labels[label.==uniqs[i]]=i
        end

        type STUD
            datas::Array{Float,2}
            baseline::Array{Float,2}
            orgsum_squares::Array{Float,3}
            orgmeans::Array{Float,2}
            orginv_cov::Array{Float,3}
            orglog_det_cov::Array{Float,1}
            ns::Array{Int,1}
            labels::Array{Int,1}
            sources::Int
            pc_max_ind::Int
            pc_gammaln_by_2::Array{Float,1}
            pc_log_pi::Float
            pc_log::Array{Float,1}
            D::Int
            N::Int
        end
         
    # precompute student-T posterior predictive distribution constants
    consts = STUD(datas,baseline,Array(Float,(D,D,max(labels))),Array(Float,(D,max(labels))),Array(Float,(D,D,max(labels))),Array(Float,max(labels)),Array(Int,max(labels)),labels,int(max(labels)),pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D,N)
