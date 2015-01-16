require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

regions = list.files('Data/HumanRegionExpr/',full.names = T)

softFile = read.design('Data/hugeHumanSoft.tsv')
regionExp = read.exp(regions[1])
corrs = matrix(data=NA, nrow = nrow(regionExp), ncol = len(regions))
colnames(corrs)=basename(regions)
rownames(corrs)=regionExp$Gene_Symbol
for (i in 1:len(regions)){
    regionExp = read.exp(regions[i])
    
    #medExpr = median(unlist(regionExp[4:len(regionExp)]))
    #medVar = median(apply(regionExp[4:len(regionExp)],1,var))
    #keep = apply(regionExp[4:len(regionExp)],1,function(row){
    #    mean(row)>medExpr | var(row)>medVar
    #})
    #regionExp = regionExp[keep,]
    
    #regionExp = mostVariable(regionExp,'Gene_Symbol')
    
    names(regionExp)[4:len(regionExp)] = sub('_.*','',names(regionExp)[4:len(regionExp)] )
    relSoft = softFile[softFile$GSM %in% names(regionExp),]
    
    corrs[,i]=cor(relSoft$pH, t(regionExp[,4:len(regionExp)])) 
}

quants = quantile(corrs,c(0.05,0.95))

keep = apply(corrs,1,function(x){max(x)>quants[2] | min(x)<quants[1]})

p=heatmap.3(corrs[keep,], trace = "none",
            col = heatPalette, cexCol=1,dendrogram = 'column')

