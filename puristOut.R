puristOut = function(geneLoc){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    fileContents = lapply(paste0(geneLoc,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    for (i in 1:length(fileContents)){
        geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V3))>0.5
                                                          & as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
    }


    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    names(puristList) = names(geneList)
    return(puristList)
}