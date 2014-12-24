require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

humanBipol = function(geneLoc, bipolLoc, bipolOut){
    dir.create(bipolOut,showWarnings = F, recursive = F)
    # adapted from Rotation 3. Very specific to the bipolar data. Fix it later
    library(reshape2)
    require(ggplot2)

    load(bipolLoc)
    bpCntScz = aned_high_GSE12649
    bpCnt = aned_high_GSE5388
    bpCntSczDes = Samples_GSE12649
    bpCntDes = Samples_GSE5388
    rm(aned_high_GSE12649)
    rm(aned_high_GSE5388)
    rm(Samples_GSE12649)
    rm(Samples_GSE5388)

    source('puristOut.R')
    puristList = puristOut(geneLoc)
    commonGround = puristList
    source('homologene.R')
    humanGenes = lapply(commonGround, function(x){mouse2human(x)$humanGene})
    
    humanGroundScz = vector(mode = 'list', length = length(commonGround))
    humanGround = vector(mode = 'list', length = length(commonGround))
    names(humanGround) = names(commonGround)
    names(humanGroundScz) = names(commonGround)
    bpCntSczExpr = bpCntScz[,4:ncol(bpCntScz)]
    bpCntSczGenes = bpCntScz[,1:3]
    colnames(bpCntSczGenes)[2] ='Gene.Symbol'
    
    bpCntExpr = bpCnt[,4:ncol(bpCnt)]
    bpCntGenes = bpCnt[,1:3]
    colnames(bpCntGenes)[2] ='Gene.Symbol'
    
    
    humanGround = lapply(humanGenes, function(x){bpCntGenes$Gene.Symbol[bpCntGenes$Gene.Symbol %in% x]})
    humanGroundScz = lapply(humanGenes, function(x){bpCntSczGenes$Gene.Symbol[bpCntGenes$Gene.Symbol %in% x]})
    

    windowSize = 14

    rownames(bpCntSczExpr) = bpCntSczGenes$Gene.Symbol
    rownames(bpCntExpr) = bpCntGenes$Gene.Symbol

    usedStuff = vector(mode = 'list', length = length(commonGround))
    names(usedStuff) = names(commonGround)

    for (i in 1:length(commonGround)){

        #bpCntSczR = t(scale(t(bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ])))
        #bpCntR = t(scale(t(bpCntExpr[which(bpCntGenes$Gene.Symbol %in% humanGround[[i]]), ])))

        #bpCntSczR= apply(bpCntSczR,2, mean)
        #bpCntR = apply(bpCntR,2, mean)
        bpCntSczR = bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ]

        indivGenes = data.frame(melt(bpCntSczR),rownames(bpCntSczR))

        indivGenes$variable = sub('_((.)|(..))','',indivGenes$variable)
        names(indivGenes) = c('group', 'expression', 'gene')

        (ggplot(indivGenes, aes(x= group, y =expression)) + geom_point() + facet_wrap('gene')  + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20), axis.title.y = element_text(vjust=0.5, size=20),axis.title.x = element_text(vjust=0.5, size=0) , title = element_text(vjust=0.5, size=20)))
        ggsave(filename = paste0(bipolOut,'/GSE12649 indivExpr ',names(commonGround)[i] , '.png'))
        #probing... get variance of genes
        varScz = apply(bpCntSczR,1,var)
        meanScz = apply(bpCntSczR,1,mean)

        bpCntR = bpCntExpr[which(bpCntGenes$Gene.Symbol %in% humanGround[[i]]), ]
        indivGenes = data.frame(melt(bpCntR),rownames(bpCntR))
        indivGenes$variable = sub('_((.)|(..))','',indivGenes$variable)
        names(indivGenes) = c('group', 'expression', 'gene')

        (ggplot(indivGenes, aes(x= group, y =expression)) + geom_point() + facet_wrap('gene') + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20), axis.title.y = element_text(vjust=0.5, size=20),axis.title.x = element_text(vjust=0.5, size=0) , title = element_text(vjust=0.5, size=20)))
        ggsave(filename = paste0(bipolOut,'/GSE5388 indivExpr ',names(commonGround)[i] , '.png'))


        varCnt = apply(bpCntR, 1, var)
        meanCnt = apply(bpCntR,1,mean)

        pcaCntScz = prcomp(t(bpCntSczR), scale = T)

        pcaCntScz$rotation = pcaCntScz$rotation * ((sum(pcaCntScz$rotation[,1])<0)*(-2)+1)
        
        write.table(pcaCntScz$rotation[,1,drop=F],
                    file = paste0(bipolOut,'/GSE12649 ',names(commonGround)[i] , '.tsv'),
                    quote = F, row.names = T, col.names = F, sep='\t')

        #report negative ones
        sczNeg = rownames(pcaCntScz$rotation)[pcaCntScz$rotation[,1]<(0)]
        sczVeryNeg = rownames(pcaCntScz$rotation)[pcaCntScz$rotation[,1]<(-0.1)]
        print(i)
        print(pcaCntScz$rotation[pcaCntScz$rotation[,1]<0,1])
        print('####')

        pcaCntScz$x = t(as.matrix(t(scale(t(bpCntSczR))))) %*% as.matrix(pcaCntScz$rotation)

        pcaCnt = prcomp(t(bpCntR), scale = T)

        pcaCnt$rotation = pcaCnt$rotation * ((sum(pcaCnt$rotation[,1])<0)*(-2)+1)
        
        write.table(pcaCnt$rotation[,1,drop=F],
                    file = paste0(bipolOut,'/GSE5388 ',names(commonGround)[i] , '.tsv'),
                    quote = F, row.names = T, col.names = F, sep = '\t')
        
        #report negative ones
        sczCnt = rownames(pcaCnt$rotation)[pcaCnt$rotation[,1]<(0)]
        print(pcaCnt$rotation[pcaCnt$rotation[,1]<0,1])

        #report common negatives
        print(intersect(sczNeg,sczCnt))


        pcaCnt$x = t(as.matrix(t(scale(t(bpCntR))))) %*% as.matrix(pcaCnt$rotation)


        bpCntSczR = t(pcaCntScz$x[,1])
        bpCntR = t(pcaCnt$x[,1])

        #bpCntSczR = bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ]
        #bpCntR = bpCntExpr[bpCntGenes$Gene.Symbol %in% humanGround[[i]], ]

        #cntSczFrame = data.frame(zScore = unlist(as.data.frame(bpCntSczR)), sample = repIndiv(bpCntSczDes['disease_state',], nrow(bpCntSczR)))
        cntSczFrame = data.frame(PC1 = unlist(as.data.frame(bpCntSczR)), sample = bpCntSczDes['disease_state',])
        
        windowSize = max(abs(cntSczFrame$PC1)) + 0.5



        #collect info about the markers
        toMerge = list(varScz = (as.numeric(varScz)),
                       varCnt = (as.numeric(varCnt)),
                       sczRot = (as.numeric(pcaCntScz$rotation[,1])),
                       CntRot = (as.numeric(pcaCnt$rotation[,1])),
                       meanScz = (as.numeric(meanScz)),
                       meanCnt = (as.numeric(meanCnt)))

        toMerge1 = list(varScz,
                        varCnt)
        names(toMerge1[[1]]) = names(varScz)
        names(toMerge1[[2]]) = names(varCnt)

        toMerge2 = list(pcaCntScz$rotation[,1],
                        pcaCnt$rotation[,1])
        names(toMerge2[[1]]) = rownames(pcaCntScz$rotation)
        names(toMerge2[[2]]) = rownames(pcaCnt$rotation)

        toMerge3 = list(meanScz,
                        meanCnt)

        names(toMerge3[[1]]) = names(meanScz)
        names(toMerge3[[2]]) = names(meanCnt)

        usedStuff[[i]] = cbind(do.call(merge, c(toMerge1,by=0, all=TRUE)),
                               do.call(merge, c(toMerge2,by=0, all=TRUE)),
                               do.call(merge, c(toMerge3,by=0, all=TRUE)))

        usedStuff[[i]]= usedStuff[[i]][,-c(4,7)]
        colnames(usedStuff[[i]]) = c('gene','varScz', 'varCnt', 'sczRot', 'cntRot', 'meanScz', 'meanCnt')

        contBP=  tryCatch({
            wilcox.test(cntSczFrame[cntSczFrame$sample=='BP',1], cntSczFrame[cntSczFrame$sample=='Cont',1])
        },error = function(e){
            return('Er')
        })
        contScz = tryCatch({
            wilcox.test(cntSczFrame[cntSczFrame$sample=='SCZ',1], cntSczFrame[cntSczFrame$sample=='Cont',1])
        },error = function(e){
            return('Er')
        })
        tryCatch({
            (p = ggplot(cntSczFrame, aes(x =sample, y =PC1  ))
             + geom_violin( color="#C4C4C4", fill="#C4C4C4")
             + geom_boxplot(width=0.1,fill = 'lightblue')
             + geom_point(size = 3)
             + ggtitle(names(commonGround)[i])
             + annotate('text' , x = 1.5, y =1, label = paste0('p = ', round(contBP$p.value,digits = 5)), size = 8)
             + annotate('text' , x = 2.5, y =1, label = paste0('p = ', round(contScz$p.value,digits = 5)), size = 8)
             + scale_y_continuous(limits=c(-windowSize, windowSize),name="Relative estimate of cell type amounts")
             + theme_bw()
             + theme(axis.text.x  = element_text(vjust=windowSize/2, size=25),
                     axis.title.y = element_text(vjust=0.5, size=25),
                     axis.title.x = element_text(vjust=0.5, size=0) ,
                     title = element_text(vjust=0.5, size=25),
                     axis.text.y = element_text(size = 13))
            )

            ggsave(filename = paste0(bipolOut,'/GSE12649 ',names(commonGround)[i] , '.png'))
        }, error = function(e){

        })

        #cntFrame = data.frame(zScore = unlist(as.data.frame(bpCntR)), sample = gsub('_t','',repIndiv(bpCntDes['disease_state',], nrow(bpCntR))))
        cntFrame = data.frame(PC1 = unlist(as.data.frame(bpCntR)), sample = gsub('_t','',bpCntDes['disease_state',]))
        windowSize = max(abs(cntFrame$PC1))+0.5
        

        contBP = tryCatch({
            wilcox.test(cntFrame[cntFrame$sample=='BP',1], cntFrame[cntFrame$sample=='Cont',1])
        },error = function(e){
            return( 'Er')
        })
        tryCatch({
            (p = ggplot(cntFrame, aes(x =sample, y =PC1  ))
             + geom_violin( color="#C4C4C4", fill="#C4C4C4")
             + ggtitle( names(commonGround)[i])
             + geom_boxplot(width=0.1,fill = 'lightblue')
             + geom_point(size = 3)
             + annotate('text' , x = 1.5, y =1, label = paste0('p = ',  round(contBP$p.value,digits = 5)), size = 8)
             + scale_y_continuous(limits=c(-windowSize, windowSize),name="Relative estimate of cell type amounts")
             + theme_bw()
             + theme(axis.text.x  = element_text(vjust=windowSize/2, size=25),
                     axis.title.y = element_text(vjust=0.5, size=25),
                     axis.title.x = element_text(vjust=0.5, size=0),
                     title = element_text(vjust=0.5, size=25), 
                     axis.text.y = element_text(size = 13))
             + scale_x_discrete(limits = c('Cont','BP'))
            )

            ggsave(filename = paste0(bipolOut, '/GSE5388 ', names(commonGround)[i] ,'.png'))
        }, error = function(e){

        })
    }
}