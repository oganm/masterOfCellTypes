source('predictAmount.R')
library(ggplot2)
library(igraph)
library(corpcor)
library(biomaRt)
library(gplots)
library(parallel)

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


# loading brain region data -------
orbitalDat = read.exp('Data/DevelopData/OFC_full')
orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
#orbitalDat = read.exp('Data/DevelopData/OFC')

list[orbGenes, orbExp] = sepExpr(orbitalDat)
meta = read.table('BrainDev/Metadata_GSE25219', sep='\t', stringsAsFactors=F)
relMeta = meta[match(cn(orbExp),rn(meta)),]
relMeta = relMeta[!duplicated(relMeta$brain.code),]
orbExp = orbExp[cn(orbExp) %in% rn(relMeta)]
orbGenes = orbGenes[!apply(orbExp,1,max)<6.5,]
orbExp = orbExp[!apply(orbExp,1,max)<6.5,]
orbitalDat = cbind(orbGenes,orbExp)

orbitalDat = mostVariable(orbitalDat)
list[orbGenes, orbExp] = sepExpr(orbitalDat)
rownames(orbExp) = orbGenes$Gene.Symbol

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


jesseCoexp = lapply(pairList, function(pairs){
    sapply(1:len(estimAge), function(i){
        
        diag(cor(t(orbExp[match(pairs[,1],orbGenes$NCBIids),order(relMeta$Age_years)[i:(i+10)]]),
                 t(orbExp[match(pairs[,2],orbGenes$NCBIids),order(relMeta$Age_years)[i:(i+10)]])))
    })
})



# plot mean coexpression for jesse genes
lapply(1:len(jesseCoexp), function(x){
    plot(y= apply(jesseCoexp[[x]], 2, mean),x = estimAge, main = names(jesseCoexp)[x])
})

lapply(1:len(jesseCoexp), function(x){
    plot(density(jesseCoexp[[x]]),main = names(jesseCoexp)[x])
})




# mean age per window ------
estimAge =  vector(length = nrow(relMeta) - 10 )

for(i in 1:len(estimAge)){
    estimAge[i] = mean(relMeta$Age_years[order(relMeta$Age_years)[i:(i+10)]])
}
plot(estimAge)
# coexpression changes ------------------------------
#coexpressions = vector(mode = 'list', length = nrow(relMeta) - 10)


# randTolerance = (max(orbExp)-min(orbExp)) * 0.05
# random gene selection function ------
generalCor = cor(t(adultExp))
limit = quantile(generalCor, .95)

meanExpr = apply(adultExp,1,mean)
meanExprQuant  =ecdf(meanExpr)

coexpSelection = apply(generalCor,2,function(x){sum(x>limit)})
coexpSelectQuant = ecdf(coexpSelection)

randomGene = function(whichGene, times = 500){
    if (is.character(whichGene)){
        whichGene = which(names(meanExpr) %in% whichGene)
}        
    # limit for mean expression
    list[minExp, maxExp] = quantile(meanExpr, 
                                    c(max(meanExprQuant(meanExpr[whichGene])-0.1,0),
                                      min(meanExprQuant(meanExpr[whichGene])+0.1,1)))
    # limit for coexpressed gene count
    list[minCoexp, maxCoexp] = quantile(coexpSelection,
                                        c(max(coexpSelectQuant(coexpSelection[whichGene])-0.1,0),
                                          min(coexpSelectQuant(coexpSelection[whichGene]) + 0.1, 1)))
    candidates = meanExpr>minExp & meanExpr<maxExp & coexpSelection>=minCoexp & coexpSelection <= maxCoexp
    whichCand = sample(which(candidates),size=times,replace=T)
    return(whichCand)
        
}


# random genes for gene lists --------
# outputs indexes
randomGenes = lapply(geneList,function(x){
    relGenes = rownames(coexpression) %in% mouse2human(x)$humanGene
    whichGenes = which(relGenes)
    sapply(whichGenes,randomGene)
})

# random genes for jesse's pairs ------------
# outputs indexes
randomPairs = mclapply(pairList,mc.cores=2,function(x){
    whichPairs = matrix(match(unlist(x),orbGenes$NCBIids), ncol = 2)
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

for(i in 1:len(estimAge)){
    print(i)
    coexpression =  cor(t(orbExp[,order(relMeta$Age_years)[i:(i+10)]]))
    limits = quantile(sm2vec(coexpression),c(seq(0,1,0.01)))
    write.table(as.data.frame(limits),paste0("Documents/Development/CoexpDevData/",round(estimAge[i]*1000)/1000), col.names=F )
    # writing data for geneLists
    for (j in 1:len(geneList)){
        dir.create(paste0('Documents/Development/CoexpDevData/GeneLists/',names(geneList)[j]), showWarnings=F)
        relGenes = rownames(coexpression) %in% mouse2human(geneList[[j]])$humanGene
        whichGenes = which(relGenes)
        relCoexp = coexpression[whichGenes,whichGenes]
        write.csv(relCoexp,paste0('Documents/Development/CoexpDevData/GeneLists/',names(geneList)[j],'/',round(estimAge[i]*1000)/1000),row.names=F)
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
                             names(randomGenes)[j], '/', round(estimAge[i]*1000)/1000),
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
                  paste0('Documents/Development/CoexpDevData/GenePairs/',names(pairList)[j], '/', round(estimAge[i]*1000)/1000),
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
                      paste0('Documents/Development/CoexpDevData/GenePairs/RandomGenes/',x,'/',names(pairList)[j], '/', round(estimAge[i]*1000)/1000),
                      row.names=F)
        })
    }
}


for (i in 1:len(estimAge)){
    
}
    

# obsolete delete later -----------------
dir.create('Documents/Development/CoexpDevelop', recursive=T)
geneCoexp = vector(mode = 'list', length = len(geneList))
 
for (i in 1:len(estimAge)){
    coexpression = cor(t(orbExp[,order(relMeta$Age_years)[i:(i+10)]] ))
    #write.csv(round(coexpression*100000)/100000,file = paste0('Documents/Development/CoexpDevData/', round(estimAge[i]*1000)/1000, '.csv'), row.names=F)
    
    limits = quantile(sm2vec(coexpression),c(0.99,0.95,0.90,1))
    print(i)
    for (j in 1:len(geneList)){
        dir.create(paste0('Documents/Development/CoexpDevelop/',names(geneList)[j]), showWarnings=F)
        relGenes = rownames(coexpression) %in% mouse2human(geneList[[j]])$humanGene
        relCoexp = coexpression[relGenes,relGenes]
        connections = relCoexp > limits['95%']
        connectionList = graph.adjacency(connections,diag=F, mode = 'undirected')
        exprLevel = apply(orbExp[relGenes,order(relMeta$Age_years)[i:(i+10)]],1,mean)
        my.colors<-colorRampPalette(c("white", "red")) #creates a function my.colors which interpolates n colors between blue, white and red
        color.df<-data.frame(COLOR_VALUE=seq(min(exprLevel),max(exprLevel),0.1),
                             color.name=my.colors(len(seq(min(exprLevel),max(exprLevel),0.1))))
        colors = color.df$color.name[
            sapply(exprLevel,function(x){
                which(abs(color.df$COLOR_VALUE-x) == min(abs(color.df$COLOR_VALUE-x)))
                })]
        V(connectionList)$color = as.character(colors)
        #V(connectionList)$label.cex = 0.5
        svg(paste0('Documents/Development/CoexpDevelop/',names(geneList)[j],'/', round(estimAge[i]*1000)/1000),
            width = 10, height = 10)
        plot(connectionList,
             #vertex.size = 14,
             layout = layout.fruchterman.reingold(connectionList, circular=T))
        legendGrad(as.char(color.df$color.name),color.df$COLOR_VALUE )
        dev.off()
        relCoexp = (relCoexp - limits['90%'])/(limits['100%'] - limits['90%'])
        diag(relCoexp) = 0
        relCoexp[relCoexp<0] = 0
        geneCoexp[[j]] =  rbind(geneCoexp[[j]] ,apply(relCoexp,2,sum))
    }
}

# mean expression in windows ------------
expressions = vector(mode = 'list', length = len(geneList))
for (i in 1:len(estimAge)){
    for (j in 1:len(geneList)){
        expressions[[j]] = cbind( expressions[[j]],
                                  apply(
                                      orbExp[rn(orbExp) %in% mouse2human(geneList[[j]])$humanGene,order(relMeta$Age_years)[i:(i+10)]],
                                      1,
                                      mean)
        )
    }
}


names(geneCoexp) = names(geneList)
# all coexp
dir.create("Documents/Development/CoexpDevelopDif")
for (i in 1:len(geneCoexp)){
    coefs  =apply(geneCoexp[[i]],2,function(x){
        toFit = data.frame(data = x, age = estimAge)
        fitted = lm(data~age, toFit)
        fitted$coefficients[2]
    })
    
    # toPlot = geneCoexp[[i]][,order(apply(geneCoexp[[i]],2,var),decreasing=T)[1:5]]
    toPlot = geneCoexp[[i]][,order(coefs,decreasing=F)]
    
    toPlot = melt(toPlot)
    toPlot$Var1 = estimAge[toPlot$Var1]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' coexpression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="Coexpression level"))
    svg(paste0('Documents/Development/CoexpDevelopDif/',names(geneCoexp[i]),'_coexp'))
    plot(p)
    dev.off()    
}



# preborn coexp
dir.create("Documents/Development/CoexpDevelopDif_PreNat")

for (i in 1:len(geneCoexp)){
    coefs  =apply(geneCoexp[[i]],2,function(x){
        toFit = data.frame(data = x[estimAge<=estimAge[18]], age = estimAge[estimAge<=estimAge[18]])
        fitted = lm(data~age, toFit)
        fitted$coefficients[2]
    })
    
    # toPlot = geneCoexp[[i]][,order(apply(geneCoexp[[i]],2,var),decreasing=T)[1:5]]
    toPlot = geneCoexp[[i]][estimAge<=estimAge[18],order(coefs,decreasing=T)]
    
    toPlot = melt(toPlot)
    toPlot$Var1 = estimAge[toPlot$Var1]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' coexpression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="Coexpression level"))
    svg(paste0('Documents/Development/CoexpDevelopDif_PreNat/',names(geneCoexp[i]),'_coexp'))
    plot(p)
    dev.off()
}

# all expressions
dir.create("Documents/Development/DevelopExpr")
for (i in 1:len(geneList)){
    toPlot = t(orbExp[rn(orbExp) %in% mouse2human(geneList[[i]])$humanGene,])
    keep = apply(toPlot,2,mean)>0
    if (sum(keep)==0){
        next
    }

    toPlot = t(orbExp[rn(orbExp) %in% mouse2human(geneList[[i]])$humanGene,])
    toPlot = toPlot[,keep,drop=F]
    toPlot = melt(toPlot)
    toPlot$Var1 = relMeta$Age_years[match(toPlot$Var1 ,rn(relMeta))]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' expression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="expression"))
    svg(paste0('Documents/Development/DevelopExpr/',names(geneCoexp[i])))
    plot(p)
    dev.off()
}

# prenatal expressions
dir.create("Documents/Development/DevelopExpr_PreNat")
for (i in 1:len(geneList)){
    toPlot = t(orbExp[rn(orbExp) %in% mouse2human(geneList[[i]])$humanGene,])
    keep = apply(toPlot,2,mean)>0
    if (sum(keep)==0){
        next
    }
    toPlot = toPlot[,keep,drop=F]
    toPlot = melt(toPlot)
    toPlot$Var1 = relMeta$Age_years[match(toPlot$Var1 ,rn(relMeta))]
    toPlot = toPlot[toPlot$Var1<estimAge[18],]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' expression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="expression"))
    svg(paste0('Documents/Development/DevelopExpr_PreNat/',names(geneCoexp[i])))
    plot(p)
    dev.off()
}
# rotation correlations ---------
estimates = vector(mode='list', length = nrow(relMeta) - 10 )
for (i in 1:len(estimates)){
    estimates[[i]]= cellTypeEstimate(exprData=cbind(orbGenes,orbExp[,order(relMeta$Age_years)[i:(i+10)] ]),
                                     genes=geneList,geneColName='Gene.Symbol',outlierSampleRemove=F)
}
dir.create('Documents/Development/PCDevelop')
for (i in 1:len(estimates[[1]]$rotations)){
    windowRot = sapply(estimates, function(x){
        x$rotations[[i]][,1]})
    windowRot = as.data.frame(windowRot)
    colnames(windowRot) = round(estimAge*1000)/1000
    tryCatch({
        heatmap.2(abs(cor(windowRot,method = 'spearman')), 
                  trace = 'none' ,
                  main = paste(names(estimates[[1]]$rotations)[[i]], 
                               nrow(estimates[[1]]$rotations[[i]])),
                  dendrogram = 'none',
                  Colv=F,
                  Rowv=F
        )
        png(paste0('Documents/Development/PCDevelop/',names(estimates[[1]]$rotations)[[i]]),width=800,height = 800)
        heatmap.2(abs(cor(windowRot,method = 'spearman')), 
                  trace = 'none' ,
                  main = paste(names(estimates[[1]]$rotations)[[i]], 
                               nrow(estimates[[1]]$rotations[[i]])),
                  dendrogram = 'none',
                  Colv=F,
                  Rowv=F
        )
        dev.off()
    },
    error = function(e){})
    #ageCor = apply(windowRot,1,cor.test,estimAge)
    #ageCor = sapply(ageCor,function(x){x$p.value})
    #toPlot = melt(windowRot[ageCor>0.05,,drop=F])
    #if (nrow(toPlot)>0){
    #    p=(ggplot(toPlot, aes(x = Var2, y = value, group = Var1)) + geom_line(aes(col = Var1)) 
    #     + ggtitle(paste(names(estimates[[1]]$rotations)[[i]],nrow(estimates[[1]]$rotations[[i]]))))
    #    plot(p)
    #}
}

# pseudo correction for mithocondrial averages
correction = apply(geneCoexp$mithocondrial,1,mean)

dir.create("Documents/Development/CoexpDevelopDif_MitoCorr")
for (i in 1:len(geneCoexp)){
    coefs  =apply(geneCoexp[[i]],2,function(x){
        toFit = data.frame(data = x, age = estimAge)
        fitted = lm(data~age, toFit)
        fitted$coefficients[2]
    })
    
    # toPlot = geneCoexp[[i]][,order(apply(geneCoexp[[i]],2,var),decreasing=T)[1:5]]
    toPlot = geneCoexp[[i]][,order(coefs,decreasing=F)]
    toPlot = toPlot / correction
    
    toPlot = melt(toPlot)
    toPlot$Var1 = estimAge[toPlot$Var1]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' coexpression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="Coexpression level"))
    svg(paste0('Documents/Development/CoexpDevelopDif_MitoCorr/',names(geneCoexp[i]),'_coexp'))
    plot(p)
    dev.off()    
}



# preborn coexp
dir.create("Documents/Development/CoexpDevelopDif_PreNat_MitoCorr")

for (i in 1:len(geneCoexp)){
    coefs  =apply(geneCoexp[[i]],2,function(x){
        toFit = data.frame(data = x[estimAge<=estimAge[18]], age = estimAge[estimAge<=estimAge[18]])
        fitted = lm(data~age, toFit)
        fitted$coefficients[2]
    })
    
    # toPlot = geneCoexp[[i]][,order(apply(geneCoexp[[i]],2,var),decreasing=T)[1:5]]
    toPlot = geneCoexp[[i]][estimAge<=estimAge[18],order(coefs,decreasing=T)]
    
    toPlot = melt(toPlot)
    toPlot$Var1 = estimAge[toPlot$Var1]
    p = (ggplot(toPlot, aes(x = Var1, y = value, color = Var2))+ geom_point()+geom_line()
         + ggtitle(paste0(names(geneCoexp)[i],' coexpression levels')) + scale_x_continuous(name='Age') + scale_y_continuous(name="Coexpression level"))
    svg(paste0('Documents/Development/CoexpDevelopDif_PreNat_MitoCorr/',names(geneCoexp[i]),'_coexp'))
    plot(p)
    dev.off()
}



# estimates of cell types, not particularly helpful ----------
estimates = cellTypeEstimate(exprData=cbind(orbGenes,orbExp),
                             genes=geneList,geneColName='Gene.Symbol',outlierSampleRemove=F)

relMeta = meta[match(names(estimates$estimates[[1]]),rn(meta)),]
# relMeta = relMeta[order(relMeta$Age_years),]

plot.new()
for(i in 1:len(estimates$estimates)){
    toPlot = data.frame(estimate = estimates$estimates[[i]], age = relMeta$Age_years)
    xrange = range(toPlot$age)
    yrange = range(toPlot$estimate)
    p = (ggplot(toPlot, aes(x= age, y = estimate))  + geom_point() + ggtitle(names(estimates$estimates)[[i]]))
    plot(p)
}
