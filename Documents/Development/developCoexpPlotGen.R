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


# mean age per window ------
estimAge =  vector(length = nrow(relMeta) - 10 )

for(i in 1:len(estimAge)){
    estimAge[i] = mean(relMeta$Age_years[order(relMeta$Age_years)[i:(i+10)]])
}


# load external connection data -------------
x=read.table(paste0('Documents/Development/CoexpDevData/AllScores/',round(estimAge[i],digits = 3)), row.names=1)
externalScore = sapply(1:len(estimAge), function(i){
    unlist(read.table(paste0('Documents/Development/CoexpDevData/AllScores/',round(estimAge[i],digits = 3)), row.names=1))
})

rownames(externalScore) = rownames(x)


# reading the distribution for gene pairs ----------------
pairScores = lapply(1:len(pairList), function(i){
    sapply(1:len(estimAge), function(j){
        data = read.table(paste0('Documents/Development/CoexpDevData/GenePairs/', names(pairList)[i],'/',round(estimAge[j],digits = 3)), header=F)
        quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
        coexpScorePair(data,quantiles)
    })
})
names(pairScores) = names(pairList)

# reading distributions for null pairs ------
pairNull =  mclapply(1:len(pairList),mc.cores=8, function(i){
    sapply(1:len(estimAge), function(j){
        quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
        mean(sapply(1:500, function(k){
            data = read.table(paste0('Documents/Development/CoexpDevData/GenePairs/RandomGenes/',k,'/', names(pairList)[i],'/',round(estimAge[j],digits = 3)), header=T)
            coexpScorePair(data,quantiles)
        }))
    })
})
names(pairNull) = names(pairList)

# plotting mean pair scores for ages, second part is for young part only --------------
# for (i in 1:len(pairList)){
#     plot(estimAge,pairScores[[i]], main = names(pairScores)[i])
#     plot(estimAge,pairNull[[i]], main = paste(names(pairScores)[i],'null'))
#     plot(estimAge, pairScores[[i]]/pairNull[[i]], main = paste(names(pairScores)[i],'normal'))
# }
# 
# for (i in 1:len(pairList)){
#     plot(estimAge[estimAge<1],pairScores[[i]][estimAge<1], main = names(pairScores)[i])
#     plot(estimAge[estimAge<1],pairNull[[i]][estimAge<1], main = paste(names(pairScores)[i],'null'))
#     plot(estimAge[estimAge<1], pairScores[[i]][estimAge<1]/pairNull[[i]][estimAge<1], main = paste(names(pairScores)[i],'normal'))
# }

# alternative normalization of paired data-------
# divides each by their own negative
pairNorm = lapply(1:len(pairList), function(i){
    sapply(1:len(estimAge), function(j){
        data = read.table(paste0('Documents/Development/CoexpDevData/GenePairs/', names(pairList)[i],'/',round(estimAge[j],digits = 3)), header=F)
        quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
        data = (data - quantiles['95%',])/(quantiles['100%',] - quantiles['95%',])
        data[data<0] = 0
        normal = sapply(1:500, function(k){
            data = read.table(paste0('Documents/Development/CoexpDevData/GenePairs/RandomGenes/',k,'/', names(pairList)[i],'/',round(estimAge[j],digits = 3)), header=T)
            quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
            data = (data - quantiles['95%',])/(quantiles['100%',] - quantiles['95%',])
            data[data<0] = 0
            data
        })
        data/rowMeans(as.data.frame(normal))
        sum(data)
    })
})
names(pairNorm) = names(pairList)
dir.create('temp')
for (i in 1:len(pairList)){
    svg(paste0('temp/', names(pairNorm)[i]))
    plot(estimAge[estimAge<1],pairNorm[[i]][estimAge<1], main = names(pairNorm)[i])
    dev.off()
}


# reading the data for geneLists -----------------
coexpScoreList = function(data,quantiles){
    data = (data - quantiles['95%',])/(quantiles['100%',] - quantiles['95%',])
    as.matrix(data) -> data
    diag(data) = 0
    data[data<0] = 0
    colSums(data)
}

listScores = lapply(1:len(geneList), function(i){
    sapply(1:len(estimAge), function(j){
        data = read.csv(paste0('Documents/Development/CoexpDevData/GeneLists/', names(geneList)[i],'/',round(estimAge[j],digits = 3)), header= T)
        quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
        coexpScoreList(data,quantiles)
    })
})
names(listScores) = names(geneList)

listNull = mclapply(1:len(geneList),mc.cores = 16, function(i){
    sapply(1:len(estimAge), function(j){
        print(paste0(i,'-stage:',j))
        quantiles = read.table(paste0('Documents/Development/CoexpDevData/',round(estimAge[j],digits = 3)), header=F,row.names=1)
        rowMeans(sapply(1:500, function(k){
            data = read.csv(paste0('Documents/Development/CoexpDevData/GeneLists/RandomGenes/',k,'/', names(geneList)[i],'/',round(estimAge[j],digits = 3)), header=T)
            coexpScoreList(data,quantiles)
        }))
    })
})
names(listNull) = names(geneList)

listNorm = lapply(1:len(listScores), function(x){listScores[[x]]/listNull[[x]]})
names(listNorm) = names(geneList)

dir.create('Documents/Development/ListNull')
for (i in 1:len(listScores)){
    normalized = melt((listNull[[i]]))
    normalized$Var2 = estimAge[match(normalized$Var2,1:len(estimAge))]
    p = ggplot(normalized, aes(x = Var2, y = value, group = Var1, color = Var1)) + 
        geom_line() + 
        ggtitle(names(listNull)[i]) + 
        scale_x_continuous('Age') + 
        scale_y_continuous('Coexpression score')
    ggsave(filename=paste0('Documents/Development/ListNull/', names(listNorm)[i],'.svg'),plot = p)
}




dir.create('Documents/Development/Raw')
for (i in 1:len(listScores)){
    normalized = melt((listScores[[i]]))
    normalized$Var2 = estimAge[match(normalized$Var2,1:len(estimAge))]
    p = ggplot(normalized, aes(x = Var2, y = value, group = Var1, color = Var1)) + 
        geom_line() + 
        ggtitle(names(listScores)[i]) + 
        scale_x_continuous('Age') + 
        scale_y_continuous('Coexpression score')
    ggsave(filename=paste0('Documents/Development/Raw/', names(listNorm)[i],'.svg'),plot = p)
}



# plotting normalized shit
dir.create('Documents/Development/Normalized')
for (i in 1:len(listNorm)){
    normalized = melt((listNorm[[i]]))
    normalized$Var2 = estimAge[match(normalized$Var2,1:len(estimAge))]
    p = ggplot(normalized, aes(x = Var2, y = value, group = Var1, color = Var1)) + 
        geom_line() + 
        ggtitle(names(listNorm)[i]) + 
        scale_x_continuous('Age') + 
        scale_y_continuous('Coexpression score')
    ggsave(filename=paste0('Documents/Development/Normalized/', names(listNorm)[i],'.svg'),plot = p)
}


# re-normalization according to pairs V1
for (i in 1:len(pairScores)){
    dir.create(paste0('Documents/Development/Normalized_',names(pairScores)[i],'_1'))
    for (j in 1:len(listNorm)){
        normalized = melt(listNorm[[j]] /(pairNorm[[i]] ))
        normalized$Var2 = estimAge[match(normalized$Var2,1:len(estimAge))]
        p = ggplot(normalized, aes(x = Var2, y = value, group = Var1, color = Var1)) + 
            geom_line() + 
            ggtitle(names(listNorm)[j]) + 
            scale_x_continuous('Age') + 
            scale_y_continuous('Coexpression score')
        ggsave(filename=paste0('Documents/Development/Normalized_',names(pairScores)[i],'_1/', names(listNorm)[j],'.svg'),plot = p)
    }
}

# re-normalization according to pairs V2
for (i in 1:len(pairScores)){
    dir.create(paste0('Documents/Development/Normalized_',names(pairScores)[i],'_2'))
    for (j in 1:len(listNorm)){
        normalized = melt(listNorm[[j]] /(pairScores[[i]]/pairNull[[i]]) )
        normalized$Var2 = estimAge[match(normalized$Var2,1:len(estimAge))]
        p = ggplot(normalized, aes(x = Var2, y = value, group = Var1, color = Var1)) + 
            geom_line() + 
            ggtitle(names(listNorm)[j]) + 
            scale_x_continuous('Age') + 
            scale_y_continuous('Coexpression score')
        ggsave(filename=paste0('Documents/Development/Normalized_',names(pairScores)[i],'_2/', names(listNorm)[j],'.svg'),plot = p)
    }
}

