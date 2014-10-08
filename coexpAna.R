require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')
require(Matrix)
require(glmnet)
coexpAna = function(humanDesLoc, genesOut, regionMapping, humanDat, coexpHOut){

    regionMap =  read.table(regionMapping,header=T,sep='\t',stringsAsFactors=F)
    allGenLocs = list.dirs(genesOut)
    humanDes = read.table(humanDesLoc,header=T,sep='\t',stringsAsFactors=F)
    # removes unnecessary studies. must match the condition used to generate files.
    humanDes = humanDes[humanDes$Region %in% regionMap[,2]
                        & humanDes$platform=='GPL5175'
                        &humanDes$Sex == 'Male',]
    humanExprs = list.files(humanDat, full.names = T)

    humanPres = vector(mode = 'list',length= len(humanExprs))
    for (i in 1:len(humanExprs)){
        humanPres[[i]] = read.csv(humanExprs[grepl(regionMap[i,2],basename(humanExprs))],
                                  header = T ,
                                  stringsAsFactors = F)
    }

    humanExpr=data.frame(1:nrow(humanPres[[i]]))
    for (i in 1:len(humanExprs)){
        humanExpr = cbind(humanExpr, humanPres[[i]][,4:len(humanPres[[i]])])
    }
    humanExpr = humanExpr[,-1]
    # removes an unnecessary part of the collumn names.
    names(humanExpr) =
        unlist(regmatches(names(humanExpr), gregexpr('GSM.*?(?=_)', names(humanExpr),perl=T)))

    humanGene = humanPres[[1]][,1:3]
    humanExpr = humanExpr[!is.na(humanGene$Gene_Symbol),]
    humanGene = humanGene[!is.na(humanGene$Gene_Symbol),]

    # filtering as described in http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0003911#s3
    # slightly modified to use variance instead of range but the spirit is the same
    medExp = median(unlist(humanExpr))
    medVar = median(apply(humanExpr,1,var))
    keep = apply(humanExpr,1,function(row){
        mean(row)>medExp | var(row)>medVar
    })
    humanExpr = humanExpr[keep,]
    humanGene = humanGene[keep,]

    print('read expression datas and did filtering')
    system('beep')
    # orthology information. this one didnt have a file. you can also use inparanoid or MetaPhOrs but don't think much will change
    source('humanMouseOrthologue.R')
    fullOrtho = humanMouseOrthologue(humanGene$Gene_Symbol)

    source('puristOut.R')

    for (i in 1:len(regionMap[,1])){

        regionExpr = humanExpr[,names(humanExpr) %in% humanDes$GSM[humanDes$Region == regionMap[i,2]]]
        # taking coexpression of everything, might be unnecessary, consider limiting region
        regionMouseLocs = allGenLocs[grepl(regionMap[i, 1], allGenLocs)]
        # uses puristOut for every region
        geneLists = lapply(regionMouseLocs, puristOut)
        names(geneLists) = basename(regionMouseLocs)

        for (j in 1:len(geneLists)){
            dir.create(paste0(coexpHOut,'/',
                       names(geneLists[j])),
                       recursive =T, showWarnings = F)
            for (k in 1:len(geneLists[[j]])){
                commons = which(humanGene$Gene_Symbol %in%
                                    fullOrtho$hgnc_symbol[fullOrtho$mgi_symbol %in% unlist(geneLists[[j]][[k]])])
                cellTypeRegionCor = cor(t(regionExpr), t(regionExpr[commons,]))
                markerCor = cellTypeRegionCor[commons,]


                markerExpr = regionExpr[commons,]
                # deal with it later. see
                #                 markerPC = prcomp(t(markerExpr), scale = T)
                #                 rots = markerPC$rotation[,1]
                #                 rots = rots * ((sum(rots)<0)*-2+1)
                #                 rots[rots>0]
                #                 plot(rots[rots>0],rep(1,len(rots[rots>0])))
                #                 # we want the ones with positive
                #                 toBundle



                tresh = quantile(unlist(cellTypeRegionCor),probs = c(0.99))
                rownames(cellTypeRegionCor) = humanGene$Gene_Symbol
                colnames(cellTypeRegionCor) = humanGene$Gene_Symbol[
                    which(humanGene$Gene_Symbol %in%
                              fullOrtho$hgnc_symbol[fullOrtho$mgi_symbol %in% unlist(geneLists[[j]][[k]])])]

                connections = vector(mode = 'list', length = ncol(cellTypeRegionCor))
                names(connections) = colnames(cellTypeRegionCor)
                for (l in 1:ncol(cellTypeRegionCor)){
                    connections[[l]] = rownames(cellTypeRegionCor)[cellTypeRegionCor[,l]>tresh]
                }

                interConnect = sapply(names(connections),findInList,connections)
                if (max(sapply(interConnect,len))){

                }

                humanList = unique(unlist(connections[interConnect[[which.max(sapply(interConnect,len))]]]))
                write.table(humanList ,
                            file =paste0(coexpHOut,'/', names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                            quote= F, col.names = F, row.names = F)

            }
        }

        # hat the fk is goingwon right ncw. ser ously R woat the fuci. are hou insane kr something?
        # are you ok now? good. I'll keep this just to remind you how hard you failed
    }
}