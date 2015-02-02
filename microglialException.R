microglialException = function(restDir,cores = 1){
    #applies to all files inside a directory recursively
    require(foreach)
    require(doMC)
    require(parallel)
    # so that I wont fry my laptop
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    
    fileNames = list.files(restDir, recursive =T )
    microGenes =  unlist(read.table('Data/microgliaList.txt'))
    fileNames = fileNames[grepl('Microglia',fileNames)]
    foreach (i = fileNames) %dopar% {
        micro = read.table(paste0(restDir,'/',i))
        micro = micro[micro$V1 %in% microGenes,]
        write.table(micro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
    }
}