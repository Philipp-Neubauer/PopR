  baseline =dlmread("baseline.csv",",",Float64)
        global baseline=baseline[2:end,2:end]'

        label = dlmread("labels.csv",",",Float64)
        label = label[2:end,2]
        label=int(label)  
        
        uniqs=unique(label,true)
        labels=Array(Int64,length(label))
        for i=1:length(uniqs)
            labels[label.==uniqs[i]]=i
        end

        type STUD
            datas::Array{Float64,2}
            baseline::Array{Float64,2}
            orgsum_squares::Array{Float64,3}
            orgmeans::Array{Float64,2}
            orginv_cov::Array{Float64,3}
            orglog_det_cov::Array{Float64,1}
            ns::Array{Int64,1}
            labels::Array{Int64,1}
            sources::Int64
            pc_max_ind::Int64
            pc_gammaln_by_2::Array{Float64,1}
            pc_log_pi::Float64
            pc_log::Array{Float64,1}
            D::Int64
            N::Int64
        end
         
    # precompute student-T posterior predictive distribution constants
    consts = STUD(datas,baseline,Array(Float64,(D,D,max(labels))),Array(Float64,(D,max(labels))),Array(Float64,(D,D,max(labels))),Array(Float64,max(labels)),Array(Int64,max(labels)),labels,max(labels),pc_max_ind,lgamma((1:pc_max_ind)/2),log(pi),log(1:pc_max_ind),D,N)
