require(RCurl)
require(ggplot2)
require(GGally)
require(lattice)
source('homologene.R')
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
# sets global variables to plot
loadCellTypes = function(correlate=F){
    allDataPre = read.csv(paste0(outFolder,'/',finalExp), header = T)
    mouseGene = allDataPre[,1:3]
    mouseExpr = allDataPre[,4:ncol(allDataPre)]
    rm(allDataPre)
    mouseDes = read.table('Data/meltedDesign.tsv',header=T,sep='\t', stringsAsFactors=F)
    mouseDes = mouseDes[match(colnames(mouseExpr),make.names(mouseDes$sampleName),),]
    rownames(mouseExpr) = mouseGene$Gene.Symbol
    
    mouseExpr <<- mouseExpr
    mouseGene <<-mouseGene
    mouseDes <<- mouseDes
    
    humanDataPre = read.csv('Data/humanBootstrap/1',header=T)
    humanGene = humanDataPre[,1:3]
    humanExpr = humanDataPre[,4:ncol(humanDataPre)]
    humanExpr = humanExpr[!is.na(humanGene$Gene_Symbol),]
    humanGene = humanGene[!is.na(humanGene$Gene_Symbol),]
    names(humanExpr) = gsub('[_].*?$','',names(humanExpr))
    humanDes = read.table('Data/hugeHumanSoft.tsv',header=T, sep='\t', stringsAsFactors=F)
    
    # humanExpr<<-humanExpr[,!colnames(humanExpr) %in% humanDes[humanDes$Region=='white matter',]$GSM]
    
    humanDes = humanDes[humanDes$GSM %in% names(humanExpr),]
    humanDes = humanDes[match(colnames(humanExpr), humanDes$GSM),] 
    
    #medExp = median(unlist(humanExpr))
    #keep = apply(humanExpr,1,max)>medExp
    #humanExpr <<- humanExpr[keep,]
    #humanGene <<- humanGene[keep,]
    
    dupGenes = names(table(humanGene$Gene_Symbol))[table(humanGene$Gene_Symbol)>1]
    
    dupExpr = humanExpr[humanGene$Gene_Symbol %in% dupGenes,]
    dupGenes = humanGene[humanGene$Gene_Symbol %in% dupGenes,]
    
    newhumanExpr = foreach (i = 1:length(unique(dupGenes$Gene_Symbol)), .combine=rbind) %dopar% {
        indexes = which(dupGenes$Gene_Symbol %in% unique(dupGenes$Gene_Symbol)[i])
        groupData = dupExpr[indexes,]
        chosen = which.max(apply(groupData,1,var))
        as.double(groupData[chosen,])
    }
    newhumanGene = dupGenes[match(unique(dupGenes$Gene_Symbol),dupGenes$Gene_Symbol),]
    colnames(newhumanExpr) = names(humanExpr)
    humanExpr =  humanExpr[!humanGene$Gene_Symbol %in% dupGenes$Gene_Symbol,]
    humanGene = humanGene[!humanGene$Gene_Symbol %in% dupGenes$Gene_Symbol,]
    humanExpr = rbind(humanExpr,newhumanExpr)
    humanGene = rbind(humanGene,newhumanGene)
    rownames(humanExpr) = humanGene$Gene_Symbol
    humanGene <<- humanGene
    humanExpr <<- humanExpr
    if (correlate){
        medExp = median(unlist(humanExpr))
        keep = apply(humanExpr,1,max)>medExp
        humanExpr = humanExpr[keep,]
        humanGene = humanGene[keep,]
        newHumanExpr = t(scale(t(humanExpr)))
        humanCor = cor(t(newHumanExpr),method='pearson')
        
        corvecs = sm2vec(humanCor)
        
        tresh = quantile(corvecs[corvecs>=0],0.99)
        humanCorBin = humanCor>treshs
        rownames(humanCorBin) = humanGene$Gene_Symbol
        colnames(humanCorBin) = humanGene$Gene_Symbol
        rownames(humanCor) = humanGene$Gene_Symbol
        colnames(humanCor) = humanGene$Gene_Symbol
        humanCorBin <<- humanCorBin
        humanCor <<- humanCor
    }
}

# give genes get plots
plotAll = function(genes, coloring, prop){
    mouseExpr = mouseExpr[,!is.na(mouseDes[,prop])]
    mouseDes = mouseDes[!is.na(mouseDes[,prop]),]
    data = mouseExpr[mouseGene$Gene.Symbol %in% genes,]
    colors = toColor(mouseDes[,prop],coloring)
    
    pairs(as.data.frame(t(data)),col = colors$cols,pch =15)

}



plotSingle = function(gene, prop, coloring, region){
    mouseExpr = mouseExpr[,!is.na(mouseDes[,prop])]
    mouseDes = mouseDes[!is.na(mouseDes[,prop]),]
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes$CellType))
    frame = data.frame(t(mouseExpr[mouseGene$Gene.Symbol %in% gene,]),mouseDes[,prop],mouseDes[,region],mouseDes$MajorType)
    if (!is.null(region)){
        names(frame) = c('gene','prop','region','Type')
    }
    if (is.null(region)){
        names(frame) = c('gene','prop','Type')
    }
    colors = toColor(mouseDes[,prop],coloring)
    manualColor = scale_colour_manual(name='prop', values = colors$palette)
    pal = colors$palette[order(names(colors$palette))]
    
    p =  ggplot(frame, aes(x = prop,y=gene))
    
    if (is.null(region)){
        p = p + geom_point(aes(color = prop),size =4)
    } else{
        p = p + geom_point(aes(color = prop,shape=region),size =4)
    }
    
    
    p = p + manualColor +
        theme(panel.background = element_rect(fill = "gray60"),legend.key = element_rect(fill = "gray60"))+
        theme(panel.grid.major=element_blank())+
        xlab(prop)+
        ylab(paste(gene,"log2 expression"))+
        theme(legend.title=element_blank()) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        theme(axis.title.x = element_text(size=20),
              axis.text.x  = element_text(angle=90, vjust=0.5, size=16)) + 
        theme(axis.title.y = element_text(size=20)) + 
        xlim(isNeuron[order(isNeuron[,1],isNeuron[,2]),2])
    
    
    if (!is.null(region) & len(unique(frame$region))>7){
        p = p +  scale_shape_manual(values=c(0:18))
        p = p +  theme(legend.box = "horizontal")
    }
    
    (p)
    
}


plotHuman = function(gene,prop='Region',mouse=F){
    if (mouse){
        gene = mouse2human(gene)$humanGene
    }
    frame = data.frame(t(humanExpr[humanGene$Gene_Symbol %in% gene,]),humanDes[,prop])
    names(frame) = c('gene','prop')
    (ggplot(frame, aes(x = prop,y=gene)) + 
         geom_point()+
         xlab(prop)+
         ylab(gene))
    
}

plotHumans = function(genes,prop='Region',mouse=F){
    if (mouse){
        genes = as.character(mouse2human(genes)$humanGene)
    }
    data = humanExpr[humanGene$Gene_Symbol %in% genes,]
    colors = toColor(humanDes[,prop])
    #split.screen(c(1,2))
   # plot.new()
   # screen(1) 
   png()
    pairs(as.data.frame(t(data)),col = colors$cols,pch =15)
    par(xpd=TRUE)
  #  screen(2); 
    legend(0,1,names(colors$palette),fill = colors$palette)
  dev.off()

}




