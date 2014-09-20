contamination = function(desFile, exprLoc, defMarkers, outDes){
    design = read.table(designLoc,header=T,sep='\t')
    allDataPre = read.csv(exprLoc, header = T)
    geneData = allDataPre[,1:3]
    exprData = allDataPre[,4:ncol(allDataPre)]
    exprData = exprData
    design = design[match(colnames(exprData),make.names(design$sampleName),),]
    defMarkers = read.table(defMarkers,header=T,sep='\t')
    defMarkers = lapply(as.list(defMarkers),trimElement,'')

    cindexes = vector(mode = 'list', length= len(defMarkers))
    names(cindexes) = names(defMarkers)


    for (i in 1:len(defMarkers)){
        newExprData = t(scale(t(exprData)))
        mi = apply(exprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),!design[,names(cindexes)[i]]],1,min)
        ma = apply(exprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),!design[,names(cindexes)[i]]],1,max)

        contaminations = apply((exprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),]-mi)/(ma-mi),2,mean)
        cindexes[[i]] = contaminations

    }
    cindexes = as.data.frame(cindexes)

    newDesign = cbind(design,cindexes)
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

    write.table(newDesign, row.names=FALSE,sep = '\t', quote=F ,file = outDes)

}

