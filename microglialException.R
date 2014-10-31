microglialException = function(restDir){
    #applies to all files inside a directory recursively
    fileNames = list.files(restDir, recursive =T )
    microGenes =  unlist(read.table('Data/microgliaList.txt'))
    fileNames = fileNames[grepl('Microglia',fileNames)]
    for (i in fileNames){
        micro = read.table(paste0(restDir,'/',i))
        micro = micro[micro$V1 %in% microGenes,]
        write.table(micro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
    }
}