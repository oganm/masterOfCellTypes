require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))

#humanMixedRegionBootstrap = function(humanSoftLoc,hCelDir, hBootOut)

#humanSoft  = read.table('Data/hugeHumanSoft.tsv', sep = '\t', header = T,quote = '')
hBootOut = 'Data/humanBootstrap'
hCelDir = 'humanRegionCel'

humanRegions = read.table('Data/hugeHumanSoft.tsv', header = T, sep = '\t', stringsAsFactors = F)

humanRegions = humanRegions[humanRegions$platform=='GPL5175'
                            & !humanRegions$Disease=='Cancer',]


require(foreach)
require(doMC)
require(parallel)
cores = 4
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)


library(oligo)
library(pd.huex.1.0.st.v2)

cels = list.celfiles(hCelDir)

dir.create(hBootOut,showWarnings = F, recursive = T)

size = 15
# foreach(i = 2:100) %dopar% {
for (i in 1:500){
    allRegions = unique(humanRegions$Region)
    subsetInd = vector(length = len(allRegions)*size)
    for (j in 1:len(allRegions)){
        subsetInd[(j*size+1-size):(j*size)] = sample(which(humanRegions$Region %in% allRegions[j]),size)
    }
    subsetRegions = humanRegions[subsetInd,]
    relevantCels = sapply(subsetRegions$GSM, function(x){which(grepl(x,cels))})
    sampleCels = cels[relevantCels]
    affyRaw = read.celfiles(paste0(hCelDir,'/',sampleCels))
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
    write.csv(aned, paste0(hBootOut,'/',i), row.names=FALSE)
}





