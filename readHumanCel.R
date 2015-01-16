library(oligo)
library(pd.huex.1.0.st.v2)
readHumanCel = function(GSMs, fileOut, humanDir){
    cels = list.celfiles(humanDir,listGzipped=T)
    whichCels = sapply(GSMs, function(x){which(grepl(x,cels))})
    sampleCels = cels[whichCels]
    affyRaw = read.celfiles(paste0(humanDir,'/',sampleCels))
    exonTS <- rma(affyRaw, target = "core")
    rm(affyRaw)
    featureData(exonTS) <- getNetAffx(exonTS, "transcript")
    #View the features of the obtained data
    # exonTS
    #Extract the expression data
    exp_value <- get("exprs", pos=assayData(exonTS))
    
    #Extract the proset ids
    PS_id <- pData(featureData(exonTS))[, "probesetid"]
    
    #Extract transcript annotations
    PS_an <- t(sapply(pData(featureData(exonTS))[, "geneassignment"], function(x) {
        annotation <- strsplit(x, " // ")[[1]][3]
        gene_symbol <- gsub(" ", "",strsplit( x, " // ")[[1]][2])
        mir <- paste(grep("MIR", strsplit(x, " // ")[[1]], value = T), collapse = " ")
        c(annotation, gene_symbol, mir)
    }))
    
    rownames(PS_an) <- PS_id
    colnames(PS_an) <-c("Annotation", "Gene_Symbol", "mRNAs")
    
    
    #Create the expression file
    aned <- merge(PS_an, exp_value, by.x="row.names", by.y="row.names", all.x=TRUE, sort=FALSE)
    rownames(aned) <- aned[,1]
    aned <- aned[,-1]
    write.csv(aned, fileOut, row.names=FALSE)
}