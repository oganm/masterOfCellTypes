require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

# temporary. till migration to run.R ------
coexpLoc = 'Data/humanBootstrap'
genesOut = 'Data/Fold/Relax'
coexSampOut = 'Data/humanBootGenes'
#####
coexpSampleAna = function(coexpLoc, genesOut, coexpSampOut, filter = T){
    allGenLocs = list.dirs(genesOut)
    allGenLocs = allGenLocs[-1]
    humanExprs = list.files(coexpLoc, full.names = T)
    humanExpr = read.csv(humanExprs[1],
                                    header = T ,
                                    stringsAsFactors = F)
    # this will be used to trim humanExprs and itself converted to humanGene
    humanGeneF = humanExpr[,1:3]
    source('humanMouseOrthologue.R')
    fullOrtho = humanMouseOrthologue(humanGeneF$Gene_Symbol)
    source('puristOut.R')
    
    for (i in 1:len(humanExprs)){
        humanExpr=read.csv(humanExprs[i],
                                      header = T ,
                                      stringsAsFactors = F)[,4:ncol(humanExpr)]
        humanExpr = humanExpr[!is.na(humanGeneF$Gene_Symbol),]
        humanGene = humanGeneF[!is.na(humanGeneF$Gene_Symbol),]

        # filtering as described in http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0003911#s3
        # slightly modified to use variance instead of range but the spirit is the same
        if (filter == T){
            medExp = median(unlist(humanExpr))
            medVar = median(apply(humanExpr,1,var))
            keep = apply(humanExpr,1,function(row){
                mean(row)>medExp | var(row)>medVar
            })
            humanExpr = humanExpr[keep,]
            humanGene = humanGene[keep,]
        }
        geneLists = lapply(allGenLocs, puristOut)
        names(geneLists) = basename(allGenLocs)
        for (j in 1:len(geneLists)){
            dir.create(paste0(coexpHOut,'/',
                              names(geneLists[j])),
                              recursive =T, showWarnings = F)
            
            
            
    }



}