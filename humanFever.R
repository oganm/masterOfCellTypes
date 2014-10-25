source('heatmap.3.r')
source('puristOut.R')
require(gplots)
require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
humanFever = function(geneLoc, bipolLoc, outLoc){
    load(bipolLoc)
    bpCntScz = aned_high_GSE12649
    bpCnt = aned_high_GSE5388
    bpCntSczDes = Samples_GSE12649
    bpCntDes = Samples_GSE5388
    rm(aned_high_GSE12649)
    rm(aned_high_GSE5388)
    rm(Samples_GSE12649)
    rm(Samples_GSE5388)
    genes = puristOut(geneLoc)
    humanGenes = lapply(genes,mouse2human)
    
    bpCntSczExpr = bpCntScz[,4:ncol(bpCntScz)]
    bpCntSczGenes = bpCntScz[,1:3]
    colnames(bpCntSczGenes)[2] ='Gene.Symbol'
    
    bpCntExpr = bpCnt[,4:ncol(bpCnt)]
    bpCntGenes = bpCnt[,1:3]
    colnames(bpCntGenes)[2] ='Gene.Symbol'
    
    geneDatas = list(GSE5388 = bpCntGenes,GSE12649 = bpCntSczGenes)
    exprDatas = list(GSE5388 = bpCntExpr, GSE12649 = bpCntSczExpr)
    dir.create(outLoc,showWarnings=F,recursive=T)
    for (i in 1:len(geneDatas)){
        for(j in 1:len(humanGenes)){
            BP = grepl('BP',colnames(exprDatas[[i]]))
            CNT = grepl('Cont',colnames(exprDatas[[i]]))
            SCZ = grepl('SCZ',colnames(exprDatas[[i]]))
            
            colors = rep('white',len(colnames(exprDatas[[i]])))
            colors[BP]='red'
            colors[CNT] = 'black'
            colors[SCZ] ='blue'
            
            png(filename = paste0(outLoc,'/',names(geneDatas)[i],'_',names(humanGenes)[j],'.png'), width = 1000, height = 1000)
            heat=heatmap.3(as.matrix(t(scale(t(exprDatas[[i]][geneDatas[[i]]$Gene.Symbol %in% humanGenes[[j]]$humanGene,])))),cexCol=1,ColSideColors=as.matrix(colors))    
            dev.off()
        }
    }
}