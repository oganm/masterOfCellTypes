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


    genes = vector()
    for (i in expandedHeat){
        filenames = list.files(paste0(genesOut,'/',i),include.dirs = FALSE)
        # if a location has less than 3 cell types, disregard. microglia is not to be trusted.
        if (len(filenames)<=elimCond){next}

        fileContents = lapply(paste0(genesOut,'/',i,'/', filenames),read.table)
        geneList = vector(mode = 'list', length = length(fileContents))
        names(geneList) = filenames
        
        
        
        if (ncol(fileContents[[1]])==3){
            # this if for a combination of fold change and silhouette coefficient
            for (i in 1:length(fileContents)){
                geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V3))>0.5
                                                                  & as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
            }
        } else if (ncol(fileContents[[1]])==1){
            # this is for a mere gene list
            for (i in 1:length(fileContents)){
                geneList[[i]] = as.character(fileContents[[i]]$V1)
            }
        } else if(ncol(fileContents[[1]])==2 & confRev){
            # this is for selection of percentages from confidence output
            for (i in 1:length(fileContents)){
                geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0.95])
            }
        } else if(ncol(fileContents[[1]])==2 & !confRev){
            # this is for selection of percentages from confidence output
            for (i in 1:length(fileContents)){
                geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))<0.05])
            }
        }else {
            stop('What kind of gibberish is this')
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