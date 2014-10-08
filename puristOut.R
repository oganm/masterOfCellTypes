puristOut = function(geneLoc){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    fileContents = lapply(paste0(geneLoc,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    if (ncol(fileContents[[1]])==3){
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V3))>0.5
                                                              & as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
        }
    } else if (ncol(fileContents[[1]])==1){
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else {
        error('What kind of gibberish is this')
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    names(puristList) = names(geneList)
    return(puristList)
}