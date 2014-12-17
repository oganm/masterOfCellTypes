heatGeneOut = function(genesOut, heatGenes,elimCond, removeOrg,confRev = F){
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


    genes = list()
    for (i in expandedHeat){
        filenames = list.files(paste0(genesOut,'/',i),include.dirs = FALSE)
        # if a location has less than 3 cell types, disregard. microglia is not to be trusted.
        if (len(filenames)<=elimCond){next}
        puristList = puristOut(paste0(genesOut,'/',i))
        genes = mergeList(genes,puristList)
    }
    return(genes)
}