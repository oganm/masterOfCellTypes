mostVariable = function(whichFile,outFile){
    require(foreach)
    require(doMC)
    require(parallel)
    cores = 4
    # so that I wont fry my laptop
    if (detectCores()<cores){ cores = detectCores()}
    registerDoMC(cores)

    allDataPre = read.csv(whichFile, header = T)


    geneData = allDataPre[,1:3]
    exprData = allDataPre[,4:ncol(allDataPre)]

    rowmax = apply(exprData, 1, max)
    discludeGenes = which(rowmax<6)
    exprData = exprData[-discludeGenes,]
    geneData = geneData[-discludeGenes,]



    newData <<- foreach (i = 1:length(unique(geneData$Gene.Symbol)), .combine=rbind) %dopar% {
        indexes = which(geneData$Gene.Symbol %in% unique(geneData$Gene.Symbol)[i])
        groupData = exprData[indexes,]
        chosen = which.max(apply(groupData,1,var))
        as.double(allDataPre[chosen,])
    }

    colnames(newData) = colnames(allDataPre)
# non parallel version
#     dim(newExprData) = c(length(unique(geneData$Gene.Symbol)), ncol(exprData))
#     c=system.time({
#      prog = txtProgressBar(min = 1, max = length(unique(geneData$Gene.Symbol)),style=3)
#
#     for (i in 1:length(unique(geneData$Gene.Symbol))){
#         indexes = which(geneData$Gene.Symbol %in% unique(geneData$Gene.Symbol)[i])
#         groupData = exprData[indexes,]
#         chosen = which.max(apply(groupData,1,var))
#         newExprData[i,] = as.double(groupData[chosen,])
#         setTxtProgressBar(prog, i)
#     }
#     close(prog)
#
#     newExprData = as.data.frame(newExprData)
#     })



    rownames(newData) = NULL

    write.csv(newData, file = outFile, row.names=FALSE)
}