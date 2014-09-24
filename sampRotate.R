sampRotate = function(designLoc,exprLoc,outLoc,groupNames, regionNames){
    require(RCurl)
    eval( expr = parse( text = getURL(
        "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
        ssl.verifypeer=FALSE)))


    design = read.table(designLoc,header=T,sep='\t')

    indexes = unique(design$originalIndex)

    forwardDes = data.frame()
    elimDes = data.frame()

    for (i in 1:len(indexes)){
        miniDes = design[design$originalIndex == indexes[i],]
        remNo = ceiling(nrow(miniDes)/3)
        whichRem = sample(1:nrow(miniDes),remNo)
        forwardDes = rbind(forwardDes, miniDes[-whichRem,])
        elimDes = rbind(elimDes, miniDes[whichRem,])
    }

    dir.create(outLoc, showWarnings = F, recursive = T)

    write.table(forwardDes, quote = F, row.names = F, col.names = T, paste0(outLoc,'/rotateForward'), sep = '\t')
    write.table(elimDes, quote = F, row.names = F, col.names = T, paste0(outLoc,'/rotateElim'), sep = '\t')

    geneSelect( paste0(outLoc,'/rotateForward'),exprLoc,outLoc,groupNames,regionNames)

}

