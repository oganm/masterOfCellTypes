# this version merges cell types into a single sample then looks for a variable.
# specific to our cell type data. not ideal but oh well...
mostVariableCT = function(whichFile,outFile,selectionNaming){

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

# this function is a generic function that looks for the most variable probeset
# of a gene. unlike the previous one, it takes in objects and outputs objects 
mostVariable = function(allDataPre,genes = 'Gene.Symbol'){
    exprData = allDataPre[,4:ncol(allDataPre)]
    rowmax = apply(exprData, 1, max)
    discludeGenes = which(rowmax<6)
    allDataPre = allDataPre[-discludeGenes,]
    exprData = exprData[-discludeGenes,]
    
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    allDataPre = allDataPre[!duplicated(allDataPre$Gene.Symbol),]
    return(allDataPre)
}