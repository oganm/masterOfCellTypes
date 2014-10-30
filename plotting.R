require(RCurl)
require(ggplot2)
require(GGally)
require(lattice)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
# sets global variables to plot
loadCellTypes <<- function(){
    allDataPre <<- read.csv(paste0(outFolder,'/',finalExp), header <<- T)
    geneData <<- allDataPre[,1:3]
    exprData <<- allDataPre[,4:ncol(allDataPre)]
    rm(allDataPre)
    design <<- read.table('Data/meltedDesign.tsv',header<<-T,sep<<-'\t')
    design <<- design[match(colnames(exprData),make.names(design$sampleName),),]
    rownames(exprData) = geneData$Gene.Symbol
}


genes= c('Neurod6','Mog','Slc6a3','Rbfox3')
genes=c('Slc6a3','Mog')
coloring = heatColors$CellType[4:len(heatColors$CellType)]
prop='CellType'
# give genes get plots
plotAll = function(genes, coloring, prop){
    exprData = exprData[,!is.na(design[,prop])]
    design = design[!is.na(design[,prop]),]
    data = exprData[geneData$Gene.Symbol %in% genes,]
    colors = toColor(design[,prop],coloring)
    
    pairs(as.data.frame(t(data)),col = colors$cols)
}





