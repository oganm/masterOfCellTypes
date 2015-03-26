# houskeeping analysis genes are selected based on low variance and high expression
# levels
require(RJSONIO)
require(gplots)
require(reshape2)
source('plotting.R')
require('VennDiagram')
require(stringr)
source('runVars.R')
# load our data ----
loadCellTypes()


# load RNA seq data -------------
# rnaSeq = read.table('Data/modGSE60361LOL', sep = '\t')
# rnaGenes = rnaSeq[2:nrow(rnaSeq),1]
# rnaExp = as.matrix(rnaSeq[2:nrow(rnaSeq),3:ncol(rnaSeq)])
# class(rnaExp) = 'numeric'
# rnaCelIDs = rnaSeq[2:nrow(rnaSeq),2]
# rnaCelIDs = as.double(as.character((rnaCelIDs)))
# colnames(rnaExp) = unlist(rnaSeq[1,2:(ncol(rnaSeq)-1),drop=T])
# rownames(rnaExp) =  rnaSeq[2:nrow(rnaSeq),1]

rnaSeq = read.table('Data/expression_mRNA_17-Aug-2014.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))

rnaExpAll = rnaExp 
# remove low expressed ones
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp > 5,]
rnaCelIDs = rnaCelIDs [maximExp > 5]
rnaExpCeil = rnaExp
rnaExpCeil[rnaExpCeil>50] = 50
# load aba data ----------------
source('../DealingWithAllenBrainAtlas/mouseRegionOntology.R')
# load shreejoy's version
regionExpression = read.csv("Data/aba_gene_exp/gene_energy_mat.csv", header =F )
regionNames = read.csv('Data/aba_gene_exp/brain_region_info.csv')
# fucking whitespace in shreejoy's version
regionNames$acronym = gsub(' $','',regionNames$acronym)
geneNames = read.csv('Data/aba_gene_exp/image_series_info.csv')
regionExpression = as.matrix(regionExpression)
colnames(regionExpression) = regionNames$acronym
rownames(regionExpression) = geneNames$acronym
regionExpressionLog = log(regionExpression, base = 2)
# filter anything below 10 since it introduces meaningless low frequency variation
regionExpressionLog[regionExpressionLog<=-10] = -10
regionExpressionLog[regionExpressionLog==-Inf] = -10
dim(regionExpressionLog[apply(regionExpressionLog,1,max)>-10,])
# get children and parents from my version
regionNames = cbind(regionNames,regions[match(regionNames$acronym,regions$acronym),
                                        c('children','parent_structure_id')])

# housekeeping genes from aba data ---------------
relevant = regionExpressionLog[,regionNames$tree_depth==9 | (regionNames$tree_depth<9 &regionNames$children=='')]
varianceRegion = apply(relevant,1,var)
quant = quantile(varianceRegion,0.05)
lowVarGenes = varianceRegion<quant
expressionRegion = apply(relevant,1,mean)
quant = quantile(expressionRegion,0.85)
highExpGenes = expressionRegion>quant
keeperOtheHouse = lowVarGenes & highExpGenes
abaKeep =  rn(relevant[keeperOtheHouse,])

# housekeeping genes from out data -------------------
relevant = mouseExpr[,!is.na(mouseDes$GabaDeep)]

variance = apply(relevant,1,var)
quant = quantile(variance,0.05)
lowVarGenes = variance<quant
expression = apply(relevant,1,mean)
quant = quantile(expression,0.85)
highExpGenes = expression>quant
keeperOtheHouse = lowVarGenes & highExpGenes
ourKeep = rn(relevant[keeperOtheHouse,])

# housekeeping genes from the paper ------------------
# selected accross many regions
houseKeep = read.table('http://www.tau.ac.il/~elieis/HKG/HK_genes.txt',sep='\t',header=F)
mouseKeep = human2mouse(gsub(' ','',houseKeep[,1]))$mouseGene

bestKeep = human2mouse(c('VPS29','VCP','SNRPD3','REEP5','RAB7A', 'PSMB4',
                          'PSMB2', 'GPI','EMC7','CHMP2A','C1orf43'))$mouseGene


# intersection of the sets --------------------------
whoKeepsTheKeepers  =list(bestKeep=  as.character(bestKeep),
                          mouseKeep = as.character(mouseKeep),
                          ourKeep = as.character(ourKeep),
                          abaKeep = as.character(abaKeep))
when = venn.diagram(whoKeepsTheKeepers[-1],filename=NULL)
plot.new()
grid.draw(when)

# enrichment of the genes ---------------------
phyper(q = len(intersect(abaKeep,mouseKeep)),
       m = len(abaKeep),
       n = nrow(regionExpressionLog)- len(abaKeep),
       k = len(intersect(mouseKeep, rn(regionExpressionLog))),
       lower.tail=F)

phyper(q = len(intersect(abaKeep,ourKeep)),
       m = len(abaKeep),
       n = nrow(regionExpressionLog)- len(abaKeep),
       k = len(intersect(ourKeep, rn(regionExpressionLog))),
       lower.tail=F)


phyper(q = len(intersect(abaKeep,ourKeep)),
       m = len(ourKeep),
       n = nrow(mouseExpr)- len(ourKeep),
       k = len(intersect(abaKeep, rn(mouseExpr))),
       lower.tail=F)


# check batch effects in RNA seq data -------------------------
genes = list(all=rn(rnaExpCeil),abaKeep=abaKeep,ourKeep=ourKeep,mouseKeep= mouseKeep)

batches = as.numeric(factor(str_extract(rnaMeta$cell_id,perl('.*?(?=_)'))))
rnaMeta$batches = batches
# batches  = as.numeric(factor(str_extract(colnames(rnaExp),perl('.*?(?=_)'))))
cors = matrix(NA,nrow = 4 , ncol = 5)
rownames(cors) = names(genes)

for (j in 1:len(genes)){
    relevant = rnaExpCeil[rownames(rnaExpCeil) %in% genes[[j]],rnaMeta$tissue %in% 'ca1hippocampus']
    keep = apply(relevant,1,function(x){max(x)>0})
    relevant= relevant[keep,]
    pcaHouse = prcomp (
        t(relevant),
        scale = T)
    
    for (i in 1:5){
        toAov = data.frame(comp = pcaHouse$x[,i], rnaMeta$batches[rnaMeta$tissue %in% 'ca1hippocampus',])
        corRes = cor.test(toAov[,1],as.numeric(factor(toAov[,2])),method='spearman')
        print(i)
        print(corRes$estimate)
        print(corRes$p.value)
        cors[j,i] = corRes$estimate
        #a=aov(comp~batches,toAov)
        #print(summary(a))
    }
}
svg('Documents/Single Cell QC/PC-MouseKeepHippocampus.svg')
heatmap.2(abs(cors[,]),trace='none',dendrogram='none',Rowv=F,Colv=F,cexRow=1,
          main = 'PC correlations to housekeeping gene lists')
dev.off()

cors = matrix(NA,nrow = 4 , ncol = 5)
rownames(cors) = names(genes)

# cortex
for (j in 1:len(genes)){
    relevant = rnaExpCeil[rownames(rnaExpCeil) %in% genes[[j]],rnaMeta$tissue %in% 'sscortex']
    keep = apply(relevant,1,function(x){max(x)>0})
    relevant= relevant[keep,]
    pcaHouse = prcomp (
        t(relevant),
        scale = T)
    
    for (i in 1:5){
        toAov = data.frame(comp = pcaHouse$x[,i], rnaMeta$batches[rnaMeta$tissue %in% 'sscortex',])
        corRes = cor.test(toAov[,1],as.numeric(factor(toAov[,2])),method='spearman')
        print(i)
        print(corRes$estimate)
        print(corRes$p.value)
        cors[j,i] = corRes$estimate
        #a=aov(comp~batches,toAov)
        #print(summary(a))
    }
}
svg('Documents/Single Cell QC/PC-MouseKeepHippocampus.svg')
heatmap.2(abs(cors[,]),trace='none',dendrogram='none',Rowv=F,Colv=F,cexRow=1,
          main = 'PC correlations to housekeeping gene lists')
dev.off()
# batch counts -----------------
counts = table(rnaMeta[rnaMeta$tissue %in% 'sscortex',c('batches','level1class')])
counts = counts[,apply(counts,2,sum)>0]
svg(paste0('Documents/Single Cell QC/','cortexAll','.svg'))
cellTypePal = toColor(cn(counts),c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
barplot(t(counts[,]),col= cellTypePal$col)
legend('topright',legend = names(cellTypePal$palette), fill = cellTypePal$palette)
dev.off()

sapply(1:ncol(counts),function(x){
    svg(paste0('Documents/Single Cell QC/','cortex-',cn(counts)[x],'.svg'))
    toPlot = rbind(counts[,x], apply(counts,1,sum)-counts[,x])
    barplot(toPlot,main =paste('cortex', cn(counts)[x]), beside = F, col = c('red','gray'),space=0.5)
    legend('topleft',legend = c('all cells',paste0()))
    dev.off()
    
    })


counts = table(rnaMeta[rnaMeta$tissue %in% 'ca1hippocampus',c('batches','level1class')])
counts=counts[,apply(counts,2,sum)>0]
svg(paste0('Documents/Single Cell QC/','hippocampus-all','.svg'))
cellTypePal = toColor(cn(counts),c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
barplot(t(counts[,]),col= cellTypePal$col)
legend('topright',legend = names(cellTypePal$palette), fill = cellTypePal$palette)
dev.off()

sapply(1:ncol(counts),function(x){
    svg(paste0('Documents/Single Cell QC/','hippocampus-',cn(counts)[x],'.svg'))
    toPlot = rbind(counts[,x], apply(counts,1,sum)-counts[,x])
    barplot(toPlot,main =paste('hippocampus', cn(counts)[x]), beside = F, col = c('red','gray'),space=0.5)
    dev.off()
    
})

