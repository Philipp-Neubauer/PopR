function elink(S)

    (n,niter)=size(S)
    
    Z=zeros(Float64,n-1,3)
    R=zeros(Float64,n+n-1,n-1)
    R[1:n,1]=1:n # R indexes leafs and groups
    MID=1:n # MID indexes current groups
    print(string("%|"))
    
    maxpair=zeros(Int,3,1)
    
    meanDI=zeros(Float64,n,n)
   
    for i=1:(n-1)
        for j=(i+1):n          
            sums=0
            for k=1:niter
                sums = sums+(S[j,k]==S[i,k])
            end
            meanDI[i,j]=meanDI[j,i]=sums/niter
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
        meanDI=[meanDI zeros(Float64,(n-s),1)]  # add a column of zeros
        meanDI=[meanDI;zeros(Float64,1,(n-s)+1)]# add a row of zeros for new group
        
        DIC=ones(n-s,1)
        
        for o=1:(n-s) # for each leaf/group calculate dist to the last merged group
            nodes=int(R[MID[o],find(R[MID[o],:].!=0)]) # nodes that make up group MID(o)
            for y=1:length(allnode)
                for z=1:length(nodes)
                    sums=0
                    for k=1:niter
                        sums = sums+int(S[allnode[y],k]==S[nodes[z],k])
                    end
                    temp = sums/niter
                    
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

        sorts=sort([ro,coll])
        meanDI = meanDI[[1:(sorts[2]-1), (sorts[2]+1):end],:]
        meanDI = meanDI[[1:(sorts[1]-1), (sorts[1]+1):end],:]
        meanDI = meanDI[:,[1:(sorts[2]-1), (sorts[2]+1):end]]
        meanDI = meanDI[:,[1:(sorts[1]-1), (sorts[1]+1):end]]
        
        MID= MID[[1:(sorts[2]-1), (sorts[2]+1):end]]
        MID= MID[[1:(sorts[1]-1), (sorts[1]+1):end]]
        
        #print(string(round(s/(n-1)*100),"% done"))
        #print(string("0%|","*"^int(round(s/(n-1)*35))," "^(35-int(round(s/(n-1)*35))),"|100% done"))
        print(string(*))
        end
        print(string("Done"))
        return(Z)

end

wd=(ARGS[1])
cd(wd)

S = int(readdlm("class_ids.csv",',')[2:end,:])

#elink(S[1:5,1:1])
Z=elink(S)

writecsv("linkages.csv",Z)
