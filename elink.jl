
if is(Int,Int64)
    typealias Float Float64
else
    typealias Float Float32
end

function mean(a::Number)
    a
end

wd=(ARGS[1])
cd(wd)

S = dlmread("class_ids.csv",",",Float)
S=int(S[2:end,2:end])

    n=size(S,1)
    Z=zeros(Float,n-1,3)
    R=zeros(Float,n+n-1,n-1)
    R[1:n,1]=1:n # R indexes leafs and groups
    MID=1:n # MID indexes current groups
    
    
    maxpair=zeros(Int,3,1)
    
    meanDI=zeros(Float,n,n)
    for i=1:(n-1)
        for j=(i+1):n
            meanDI[i,j]=mean(int(S[j,:].==S[i,:]))
            meanDI[j,i]=meanDI[i,j]
            if meanDI[i,j]>maxpair[3]
                maxpair=[i j meanDI[i,j]]
            end
        end
    end
    Z[1,:]=[maxpair[1:2], 1-maxpair[3]]
    
    R[n+1,1:2]= Z[1,1:2]

    meanDI= meanDI[[1:(int(Z[1,2])-1), (int(Z[1,2])+1):end],:]
    meanDI= meanDI[[1:(int(Z[1,1])-1), (int(Z[1,1])+1):end],:]
    
    meanDI= meanDI[:,[1:(int(Z[1,2])-1), (int(Z[1,2])+1):end]]
    meanDI= meanDI[:,[1:(int(Z[1,1])-1), (int(Z[1,1])+1):end]]
    
    MID=MID[[1:(int(Z[1,2])-1), (int(Z[1,2])+1):end]]
    MID=MID[[1:(int(Z[1,1])-1), (int(Z[1,1])+1):end]]

    
    for s = 2:(n-1)
        
        an1=int(R[int(Z[s-1,1]),find(R[int(Z[s-1,1]),:].!=0)])
        an2=int(R[int(Z[s-1,2]),find(R[int(Z[s-1,2]),:].!=0)])
        allnode = [an1 an2] #just merged nodes
        MID = [MID, n+(s-1)] #new group
        meanDI=[meanDI zeros(Float,(n-s),1)]  # add a column of zeros
        meanDI=[meanDI;zeros(Float,1,(n-s)+1)]# add a row of zeros for new group
        
        DIC=ones(n-s,1)
        
        for o=1:(n-s) # for each leaf/group calculate dist to the last merged group
            nodes=int(R[MID[o],find(R[MID[o],:].!=0)]) # nodes that make up group MID(o)
            for y=1:length(allnode)
                for z=1:length(nodes)
                    
                    temp = mean(mean(int(S[allnode[y],:].==S[nodes[z],:])))
                    
                    if temp<DIC[o]
                        DIC[o]=temp
                    end
                end
            end
        end
            
        meanDI[1:(n-s),(n-s)+1]=DIC
        meanDI[(n-s)+1,1:(n-s)]=DIC'
            
        ro=1
        coll=1
        for i=1:(n-s)
            for j=i:(n-s+1)
                if meanDI[i,j]>=meanDI[ro,coll]
                    ro=i
                    coll=j
                end
            end
        end  
            
        Z[s,:] = [MID[ro] MID[coll] 1-meanDI[ro,coll]] 
            
        newgr1=R[MID[ro],find(R[MID[ro],:].!=0)]
        newgr2=R[MID[coll],find(R[MID[coll],:].!=0)]
            
        
        if s.==(n-1)
            break
        end

        R[n+s,1:(length(newgr1)+length(newgr2))]=[newgr1 newgr2]

        sorts=sort([ro coll])
        meanDI = meanDI[[1:(sorts[2]-1), (sorts[2]+1):end],:]
        meanDI = meanDI[[1:(sorts[1]-1), (sorts[1]+1):end],:]
        meanDI = meanDI[:,[1:(sorts[2]-1), (sorts[2]+1):end]]
        meanDI = meanDI[:,[1:(sorts[1]-1), (sorts[1]+1):end]]
        
        MID= MID[[1:(sorts[2]-1), (sorts[2]+1):end]]
        MID= MID[[1:(sorts[1]-1), (sorts[1]+1):end]]
           
            
        end

csvwrite("linkages.csv",Z)
