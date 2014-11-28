mostVariable = function(whichFile,outFile,selectionNaming){

    allDataPre = read.csv(whichFile, header = T)
    design = read.design('Data/meltedDesign.tsv')


    exprData = allDataPre[,4:ncol(allDataPre)]
    
    cellTypes = trimNAs(unique(design[,selectionNaming]))
    
    cellTypeExpr = lapply(cellTypes,function(x){
        apply(exprData[,design[,selectionNaming] %in% x,drop=F],1,mean)
    })
    exprData = as.data.frame(cellTypeExpr)

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