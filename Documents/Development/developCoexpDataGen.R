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


# loading and expanding gene list -------
print(getwd())
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



# loading adult cortex data for random gene selection purposes ------
adultDat = read.exp('Data/HumanRegionExpr/frontal cortex')
list[adultGenes, adultExp] = sepExpr(adultDat)
adultExp = adultExp[match(orbGenes$Probe,adultGenes$Probe),]
adultGenes = adultGenes[match(orbGenes$Probe,adultGenes$Probe),]

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





# plot mean coexpression for jesse genes------
# 
# jesseCoexp = lapply(pairList, function(pairs){
#     sapply(1:len(estimAge), function(i){
#         
#         diag(cor(t(orbExp[match(pairs[,1],orbGenes$NCBIids),order(relMeta$Age_years)[i:(i+10)]]),
#                  t(orbExp[match(pairs[,2],orbGenes$NCBIids),order(relMeta$Age_years)[i:(i+10)]])))
#     })
# })

# lapply(1:len(jesseCoexp), function(x){
#     plot(y= apply(jesseCoexp[[x]], 2, mean),x = estimAge, main = names(jesseCoexp)[x])
# })
# 
# lapply(1:len(jesseCoexp), function(x){
#     plot(density(jesseCoexp[[x]]),main = names(jesseCoexp)[x])
# })




# mean age per window ------
estimAge =  vector(length = nrow(relMeta) - 10 )

for(i in 1:len(estimAge)){
    estimAge[i] = mean(relMeta$Age_years[order(relMeta$Age_years)[i:(i+10)]])
}
ageRange = vector(mode = 'list' , length = nrow(relMeta) - 10 )
for(i in 1:len(estimAge)){
    ageRange[[i]] = range(relMeta$Age_years[order(relMeta$Age_years)[i:(i+10)]])
}
# estimAge = estimAge[estimAge<5]
# coexpression changes ------------------------------
#coexpressions = vector(mode = 'list', length = nrow(relMeta) - 10)


# random gene selection function ------
generalCor = cor(t(adultExp))
limit = quantile(generalCor, .95)

meanExpr = apply(adultExp,1,mean)
meanExprQuant  =ecdf(meanExpr)

coexpSelection = apply(generalCor,2,function(x){sum(x>limit)})
coexpSelectQuant = ecdf(coexpSelection)


randomCandidate = function(whichGene){
        if (is.character(whichGene)){
            whichGene = which(names(meanExpr) %in% whichGene)
        } 
        list[minExp, maxExp] = quantile(meanExpr, 
                                        c(max(meanExprQuant(meanExpr[whichGene])-0.1,0),
                                          min(meanExprQuant(meanExpr[whichGene])+0.1,1)))
        # limit for coexpressed gene count
        list[minCoexp, maxCoexp] = quantile(coexpSelection,
                                            c(max(coexpSelectQuant(coexpSelection[whichGene])-0.1,0),
                                              min(coexpSelectQuant(coexpSelection[whichGene]) + 0.1, 1)))
        candidates = meanExpr>minExp & meanExpr<maxExp & coexpSelection >= minCoexp & coexpSelection <= maxCoexp
        return(candidates)
        # rep(T,len(meanExpr)) # to return all true
}



# random candidate sufficiency testing 
# plot(density(sapply(1:nrow(orbGenes),function(x){sum(randomCandidate(x))})))

randomGene = function(whichGene, times = 500){
    candidates = randomCandidate(whichGene)
    whichCand = sample(which(candidates),size=times,replace=T)
    return(whichCand)
}
# random genes for gene lists --------
# outputs indexes
randomGenes = lapply(geneList,function(x){
    relGenes = rownames(orbExp) %in% mouse2human(x)$humanGene
    whichGenes = which(relGenes)
    sapply(whichGenes,randomGene)
})

# random genes for jesse's pairs ------------
# outputs indexes
randomPairs = mclapply(1:len(pairList),mc.cores=2,function(x){
    print(paste0("randomPairs:", x)
    whichPairs = matrix(match(unlist(pairList[[x]]),orbGenes$NCBIids), ncol = 2)
    indivs = unique(as.double(whichPairs))
    randoms = sapply(indivs, randomGene)
    return(apply(randoms, 1, function(y){
        as.data.frame(apply(whichPairs,c(1,2), function(z){
            y[indivs %in% z]
        }))
    }))   
})

# getting coexpression data for windows (gene lists, pairs and random shit) ----------------
dir.create(paste0('Documents/Development/CoexpDevData/'), showWarnings = F)
dir.create(paste0('Documents/Development/CoexpDevData/GeneLists'), showWarnings = F)
dir.create(paste0('Documents/Development/CoexpDevData/GenePairs'), showWarnings = F)
dir.create(paste0('Documents/Development/CoexpDevData/AllScores'), showWarnings = F)

for(i in 1:len(estimAge)){
    print(paste0('mainProcess:',i))
    coexpression =  cor(t(orbExp[,order(relMeta$Age_years)[i:(i+10)]]))
    limits = quantile(sm2vec(coexpression),c(seq(0,1,0.01)))
    write.table(as.data.frame(limits),paste0("Documents/Development/CoexpDevData/",round(estimAge[i],digits = 3)), col.names=F )
    
    # total coexpression score
    # lim = read.table(paste0("Documents/Development/CoexpDevData/",round(estimAge[i],digits = 3)))
    # limits = lim[,2]
    # names(limits) = lim[,1]
    coexpScore = (coexpression - limits['95%'])/(limits['100%'] - limits['95%'])
    coexpScore[coexpScore<0] = 0
    diag(coexpScore) = 0
    coexpScore = apply(coexpScore,2,sum)
    write.table(as.data.frame(coexpScore),paste0("Documents/Development/CoexpDevData/AllScores/",round(estimAge[i],digits = 3)), col.names=F )
    
    # writing data for geneLists
    for (j in 1:len(geneList)){
        dir.create(paste0('Documents/Development/CoexpDevData/GeneLists/',names(geneList)[j]), showWarnings=F)
        relGenes = rownames(coexpression) %in% mouse2human(geneList[[j]])$humanGene
        whichGenes = which(relGenes)
        relCoexp = coexpression[whichGenes,whichGenes]
        write.csv(relCoexp,paste0('Documents/Development/CoexpDevData/GeneLists/',names(geneList)[j],'/',round(estimAge[i],digits = 3)),row.names=F)
    }
    # writing data for randomized genelists
    for(j in 1:len(randomGenes)){
        daRandom = randomGenes[[j]]
        sapply(1:nrow(daRandom),function(x){
            whichGenes = daRandom[x,]
            relCoexp = coexpression[whichGenes,whichGenes]
            dir.create(paste0('Documents/Development/CoexpDevData/GeneLists/RandomGenes/',x,'/',names(randomGenes)[j], '/'),
                       recursive=T,
                       showWarnings=F)
            write.csv(relCoexp,
                      paste0('Documents/Development/CoexpDevData/GeneLists/RandomGenes/',
                             x,'/',
                             names(randomGenes)[j], '/', round(estimAge[i],digits = 3)),
                      row.names=F)
        })
    }
    
    # writing data for gene pairs
    for (j in 1:len(pairList)){
        whichPairs = matrix(match(unlist(pairList[j]),orbGenes$NCBIids), ncol = 2)
        relCoexp = apply(whichPairs,1,function(x){
            coexpression[x[1],x[2]]
        })
        dir.create(paste0('Documents/Development/CoexpDevData/GenePairs/',names(pairList)[j]), recursive=T,showWarnings=F)
        write.table(relCoexp,
                    paste0('Documents/Development/CoexpDevData/GenePairs/',names(pairList)[j], '/', round(estimAge[i],digits = 3)),
                    row.names=F,col.names=F)
    }
    
    # for writing data for randomized gene pairs
    for (j in 1:len(randomPairs)){
        lapply(1:len(randomPairs[[j]]), function(x){
            relCoexp = apply(randomPairs[[j]][[x]],1,function(x){
                coexpression[x[1],x[2]]  
            })
            dir.create(paste0('Documents/Development/CoexpDevData/GenePairs/RandomGenes/',x,'/',names(pairList)[j]),
                       recursive=T,
                       showWarnings=F)
            write.csv(relCoexp,
                      paste0('Documents/Development/CoexpDevData/GenePairs/RandomGenes/',x,'/',names(pairList)[j], '/', round(estimAge[i],digits = 3)),
                      row.names=F)
        })
    }
}