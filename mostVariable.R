mostVariable = function(whichFile,outFile,cores = 16){
    require(foreach)
    require(doMC)
    require(parallel)
    # so that I wont fry my laptop
    if (detectCores()<cores){ cores = detectCores()}
    registerDoMC(cores)

    allDataPre = read.csv(whichFile, header = T)


    exprData = allDataPre[,4:ncol(allDataPre)]

    rowmax = apply(exprData, 1, max)
    discludeGenes = which(rowmax<6)
    allDataPre = allDataPre[-discludeGenes,]
    exprData = exprData[-discludeGenes,]
    
    # you bloody idiot.... taken from lila
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    allDataPre = allDataPre[!duplicated(allDataPre$Gene.Symbol),]
    


    write.csv(allDataPre, file = outFile, row.names=FALSE)
}