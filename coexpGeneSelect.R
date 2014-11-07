require(foreach)
require(doMC)
require(parallel)
require(corpcor)
cores = 10
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)

coexpTheVeryEnd = 'Data/CoexpGenes'
loopAround = list.dirs('Data/CoexpNetwork')[-1]

foreach (i = loopAround) %dopar%{
    files = list.files(i)
    dir.create(paste0(coexpTheVeryEnd,'/',basename(i)),recursive=T,showWarnings=F)
    for (j in files){
        network = read.csv(paste0(i,'/',j),header=T, row.names=1)
        
        confNetwork = network > 0.95
        connections = vector(mode = 'list', length = ncol(confNetwork))
        
        names(connections) = colnames(confNetwork)
        for (l in 1:ncol(confNetwork)){
            connections[[l]] = rownames(confNetwork)[confNetwork[,l]==1]
        }
        
        interConnect = sapply(names(connections),findInList,connections)
        if (max(sapply(interConnect,len))<=1){
            file.create(paste0(coexpTheVeryEnd,'/',basename(i),'/',j))
            next 
        }
        
        #this part is to be able to limit genes that just happen to connect to one gene. but I won't be filtering anything for now
        connectionCounts=table(unlist(connections[interConnect[[which.max(sapply(interConnect,len))]]]))
        humanList = unique(unlist(connections[interConnect[[which.max(sapply(interConnect,len))]]]))
        write.table(humanList ,
                    file =paste0(coexpTheVeryEnd,'/',basename(i),'/',j) ,
                    quote= F, col.names = F, row.names = F)
    }
}
