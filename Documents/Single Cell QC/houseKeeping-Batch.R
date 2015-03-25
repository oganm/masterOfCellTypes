# houskeeping analysis genes are selected based on low variance and high expression
# levels
require(RJSONIO)
require(gplots)
require(reshape2)
source('plotting.R')
require('VennDiagram')
# load our data ----
loadCellTypes()


# load RNA seq data -------------
rnaSeq = read.table('Data/modGSE60361LOL', sep = '\t')
rnaGenes = rnaSeq[2:nrow(rnaSeq),1]
rnaExp = as.matrix(rnaSeq[2:nrow(rnaSeq),3:ncol(rnaSeq)])
class(rnaExp) = 'numeric'
rnaCelIDs = rnaSeq[2:nrow(rnaSeq),2]
rnaCelIDs = as.double(as.character((rnaCelIDs)))
colnames(rnaExp) = unlist(rnaSeq[1,2:(ncol(rnaSeq)-1),drop=T])
rownames(rnaExp) =  rnaSeq[2:nrow(rnaSeq),1]
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



batches  = as.numeric(factor(str_extract(colnames(rnaExp),perl('.*?(?=_)'))))
cors = matrix(NA,nrow = 4 , ncol = 5)
rownames(cors) = names(genes)

for (j in 2:len(genes)){
    pcaHouse = prcomp (
        t(rnaExpCeil[rownames(rnaExpCeil) %in% genes[[j]],]),
        scale = T)
    
    for (i in 1:5){
        toAov = data.frame(comp = pcaHouse$x[,i], batches)
        corRes = cor.test(toAov[,1],as.numeric(factor(toAov[,2])),method='spearman')
        print(i)
        print(corRes$estimate)
        print(corRes$p.value)
        cors[j,i] = corRes$estimate
        #a=aov(comp~batches,toAov)
        #print(summary(a))
    }
}
svg('PC-MouseKeep')
heatmap.2(abs(cors[,]),trace='none',dendrogram='none',Rowv=F,Colv=F,cexRow=1,
          main = 'PC correlations to housekeeping gene lists')
dev.off()
