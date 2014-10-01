microglialException = function(restDir){
    groupNames = list.files(restDir)
    microGenes =  unlist(read.table('Data/microgliaList.txt'))
    for (i in groupNames){
        if (!'Microglia' %in% list.files(paste0(restDir,'/',i))){
            next
        }
        micro = read.table(paste0(restDir,'/',i,'/Microglia'))
        micro = micro[micro$V1 %in% microGenes,]
        write.table(micro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'/Microglia'))
    }
}