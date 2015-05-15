source('predictAmount.R')
source('mostVariable.R')
source('puristOut.R')
source('runVars.R')
library(ggplot2)
library(igraph)
library(corpcor)
library(biomaRt)
library(gplots)
library(parallel)
require(RCurl)
library(reshape2)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
# loading brain region data -------
orbitalDat = read.exp('Data/DevelopData/OFC_full')
orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
#orbitalDat = read.exp('Data/DevelopData/OFC')

list[orbGenes, orbExp] = sepExpr(orbitalDat)
meta = read.table('BrainDev/Metadata_GSE25219', sep='\t', stringsAsFactors=F)
relMeta = meta[match(cn(orbExp),rn(meta)),]
relMeta = relMeta[!duplicated(relMeta$brain.code),]
orbExp = orbExp[cn(orbExp) %in% rn(relMeta)]

# filtering for age
orbExp = orbExp[,relMeta$Age_years<5]
relMeta = relMeta[relMeta$Age_years<5,]

orbGenes = orbGenes[!apply(orbExp,1,max)<6.5,]
orbExp = orbExp[!apply(orbExp,1,max)<6.5,]
orbitalDat = cbind(orbGenes,orbExp)

orbitalDat = mostVariable(orbitalDat)
list[orbGenes, orbExp] = sepExpr(orbitalDat)
rownames(orbExp) = orbGenes$Gene.Symbol


coexpScorePair = function(coexpScorePair, quantiles){
    data = (data - quantiles['95%',])/(quantiles['100%',] - quantiles['95%',])
    data[data<0] = 0
    sum(data)
}



# loading and expanding gene list -------
geneList = puristOut('Data/RotSel/Relax/Cortex_PyramidalDeep/')
library(GO.db)

GOlist = as.list(GOTERM)

goterm  =sapply(GOlist, function(x){
    out = c(ID = GOID(x), term = Term(x),definition = Definition(x), ontology = Ontology(x))
    
    return(out)
})

goterm = as.data.frame(goterm)
goterm = t(goterm)

myelinationTerms = goterm[grep('myelination',goterm[,'term']),]

splitGo = strsplit(orbGenes$GOTerms,'[|]')

myelinGenes = orbGenes$Gene.Symbol[
    sapply(splitGo,function(x){
        any(x %in% myelinationTerms[,'ID'])
    })
    ]


geneList$OligoMyelin = geneList$Oligo[geneList$Oligo %in% mouse2human(geneList$Oligo)$mouseGene[mouse2human(geneList$Oligo)$humanGene %in% myelinGenes]]
geneList$Myelin = human2mouse(myelinGenes)$mouseGene
geneList$OligoNoMyelin = geneList$Oligo[!geneList$Oligo %in% mouse2human(geneList$Oligo)$mouseGene[mouse2human(geneList$Oligo)$humanGene %in% myelinGenes]]


mithocondrialTerms = goterm[grep('mitochondria',goterm[,'term']),]

mithocondrialGenes = orbGenes$Gene.Symbol[
    sapply(splitGo,function(x){
        any(x %in% mithocondrialTerms[,'ID'])
    })
    ]

geneList$mithocondrial =  human2mouse(mithocondrialGenes)$mouseGene


glycolysisTerms = goterm[grep('positive regulation of glycolysis',goterm[,'term']),,drop=F]
glycolysisGenes = orbGenes$Gene.Symbol[
    sapply(splitGo,function(x){
        any(x %in% glycolysisTerms[,'ID'])
    })
    ]



# positive controls from Jesse G. ----------------
load('../AuPairWise/data/pairs.Rdata')

pairList = list(stoich = stoich.pairs,
                ppin = ppin.pairs,
                yeastmouse.bottom = yeastmouse.bottom.pairs,
                yeastmouse.top = yeastmouse.top.pairs,
                multiple.complexes = multiple.complexes.pairs[,2:3],
                single.complexes = single.complexes.pairs[,2:3])

# filtering for genes that are in our chip
pairList = lapply(pairList, function(pairs){
    pairs[apply(pairs,1,function(x){
        all(x %in% orbGenes$NCBIids)
    }),]
})


estimates = cellTypeEstimate(exprData=cbind(orbGenes,orbExp),
                             genes=geneList,geneColName='Gene.Symbol',outlierSampleRemove=F)

lapply(1:len(estimates$estimates),function(i){
    plot(relMeta$Age_years,estimates$estimates[[i]], main = names(estimates$estimates)[i])
    abline(v  = relMeta$Age_years[relMeta$age %in% '35 PCW'], col='red')
})


lapply(1:len(geneList), function(i){
    toPlot = orbExp[rn(orbExp) %in% mouse2human(geneList[[i]])$humanGene,]
    toPlot$gene= rn(toPlot)
    toPlot = melt(toPlot)
    toPlot$Age = relMeta$Age_years[match(toPlot$variable, rn(relMeta))]
    toPlot = toPlot[toPlot$Age<=1,]
    (ggplot(toPlot, aes(x=Age, y = value)) + geom_line(aes(color = gene)) + ggtitle(names(geneList)[i]))
})

lapply(1:len(geneList), function(i){
    toPlot = orbExp[rn(orbExp) %in% mouse2human(geneList[[i]])$humanGene,]
    toPlot = apply(toPlot,2,mean)
    ages = relMeta$Age_years[match(names(toPlot), rn(relMeta))]
    plot(ages,toPlot, main =names(geneList)[i] )
})



