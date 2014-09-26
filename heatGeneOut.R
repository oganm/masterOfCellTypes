heatGeneOut = function(genesOut, heatGenes,elimCond, removeOrg){
    # selects genes based on a given naming scheme, including genes that aere
    # found based on regions, removes groups with less than or equal to cell
    # types with elimCond

    expandedHeat = vector()
    for (i in heatGenes){
        expandedHeat = c(expandedHeat,
                         list.files(genesOut,include.dirs = T) [grepl(paste0('((_)|(^))',i),list.files(genesOut,include.dirs = T))]
        )
    }

    # remove original folder, not based on location
    if (removeOrg) {expandedHeat = expandedHeat[-which(expandedHeat==heatGenes)]}


    genes = vector()
    for (i in expandedHeat){
        filenames = list.files(paste0(genesOut,'/',i),include.dirs = FALSE)
        # if a location has less than 3 cell types, disregard. microglia is not to be trusted.
        if (len(filenames)<=elimCond){next}

        fileContents = lapply(paste0(genesOut,'/',i,'/', filenames), read.table)
        geneList = vector(mode = 'list', length = length(fileContents))
        names(geneList) = filenames
        for (j in 1:length(fileContents)){
            geneList[[j]] = as.character(fileContents[[j]]$V1[(as.numeric(as.character(fileContents[[j]]$V3))>0.5)&(as.numeric(as.character(fileContents[[j]]$V2))>0)])
        }

        puristList = vector(mode = 'list', length = length(geneList))
        for (i in 1:length(geneList)){
            puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
        }
        genes = c(genes, unlist(puristList))
    }
    genes= unique(genes)
    return(genes)

}