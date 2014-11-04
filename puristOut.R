require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))

puristOut = function(geneLoc, lilah = F){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    fileContents = lapply(paste0(geneLoc,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    
    if (ncol(fileContents[[1]])==3 & lilah == F){
        # this if for a combination of fold change and silhouette coefficient
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V3))>0.5
                                                              & as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
        }
    }else if (ncol(fileContents[[1]])==3 & lilah == T){
        # this if for lilah's selection method
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V3))*
                                                              as.numeric(as.character(fileContents[[i]]$V2))>2])
        }
    } else if (ncol(fileContents[[1]])==1){
        # this is for a mere gene list
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else if(ncol(fileContents[[1]])==2){
        # this is for selection of percentages from confidence output
        for (i in 1:length(fileContents)){
        geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0.95])
        }
    }else {
        stop('What kind of gibberish is this')
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    names(puristList) = names(geneList)
    return(puristList)
}