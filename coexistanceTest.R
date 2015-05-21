# coexistance of genes
require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

# load rna seq data--------------
rnaSeq = read.table('Data/RNASeq/expression_mRNA_17-Aug-2014.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
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
rnaCelIDs = rnaCelIDs[maximExp > 5]
tresholds = read.table('Data/RNASeq/tresholds')

presenceAll = t(sapply(1:nrow(rnaExpAll),function(x){
    rnaExpAll[x,]>=tresholds[x,2]
}))
 presenceAll = rnaExpAll>0
rownames(presenceAll) = rownames(rnaExpAll)

probs = apply(presenceAll,1, function(x){sum(x)})
probs = sort(probs)

ourGenes = puristOut('Data/RotSel/Relax/Cortex_PyramidalDeep/')
# get real coexistance prevelance ------------------
realProbs = 
    sapply(ourGenes,function(x){
        if (len(x)<2){
            return(NA)
        }
        relevant = presenceAll[rn(presenceAll) %in% x,]
        relProbs = probs[names(probs) %in% x]
        relProbs = relProbs[match(rn(relevant), names(relProbs))]
        mean(apply(relevant,2,function(y){
            if (sum(y)>=len(y)/3){
             sum(y/relProbs)#*sum(y)
            } else{
                0
            }
            }))
        })



names(realProbsFine) = names(ourGenes)
ourGenesEx = lapply(ourGenes, function(x){
    x[x %in% names(probs)]
})

# simulate coexistance prevelance ------------------

simuProbs = 
    sapply(1:len(ourGenesEx),function(x){
        print(x)
        if (len(ourGenesEx[[x]])<2){
            return(NA)
        }
        simuGenes = sapply(ourGenesEx[[x]], function(y){
            prob = probs[names(probs) == y]
            loc = which(names(probs) == y)
            # eligible = (loc-500):(loc+500)
            eligible = which((probs < prob + prob*0.2) & (probs > prob - prob*0.2))
            
            return(names(probs[sample(eligible,500,replace=T)]))
        })
        
        return(apply(simuGenes,1,function(y){
            relevant = presenceAll[rn(presenceAll) %in% y,]
            relProbs = probs[names(probs) %in% y]
            relProbs = relProbs[match(rn(relevant), names(relProbs))]
            mean(apply(relevant,2,function(z){
                if (sum(z)>=len(z)/3){
                    sum(z/relProbs)#*sum(z)
                } else{
                    0
                }
            })
                     )
        }))
    })
names(simuProbs) = names(ourGenesEx)

# plot simulation results -------------------
dir.create('Data/SingleCelCoExists')

lapply(1:ncol(simuProbs),function(x){
    if (is.na(simuProbs[[x]])){
        return(NULL)
    }
    svg(paste0('Data/SingleCelCoExists/', names(simuProbs)[x],'.svg'))
    plot(density(simuProbs[,x]), main= names(simuProbs)[x],
         xlim=c(min(min(density(simuProbs[,x])$x),  realProbs[x]),max(max(density(simuProbs[,x])$x), realProbs[x])),
        xlab='',cex.axis = 1.3, cex.main = 1.5, cex.lab = 1.5)
   abline(v=realProbs[x],col='red', lwd = 3)
  dev.off()
})


# heatmaps of existance --------------
require(gplots)
for (i in 1:len(ourGenes)){
    toPlot = (matrix(as.numeric(presenceAll[rownames(presenceAll) %in% ourGenes[[i]],]),ncol=ncol(presenceAll)))
    if (nrow(toPlot)<=1){
        next
    }
    rownames(toPlot) = rownames(presenceAll)[rownames(presenceAll) %in% ourGenes[[i]]]
    png(paste0("Data/SingleCelCoExists/", names(ourGenes)[i],'_heat.png'), height = 800, width= 800)
    heatmap.2(t(toPlot),trace= 'none', Rowv=F, Colv=F,dendrogram='none',main = names(ourGenes)[i])
    dev.off()
}

