#designLoc = 'Data2/normalizedDesign'
#exprLoc = 'Data2/mostVariableQuantileNormalized'
#outLoc = 'Data2/Rest'
#collumn names that define groups
#groupNames = c('someNaming2',
#               'someNaming1',
#               'ourNaming3',
#               'ourNaming2',
#               'ourNaming1')

#groupNames = 'someNaming'

geneSelect = function(designLoc,exprLoc,outLoc,groupNames, regionNames){
    require(RCurl)
    eval( expr = parse( text = getURL(
        "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
        ssl.verifypeer=FALSE) ))
    require(reshape)

    require(cluster)

    #gene selector, outputs selected genes and their fold changes
    foldChange = function (group1, group2, f = 10){


        groupAverage1 = group1



        groupAverage2 = tryCatch({apply(group2, 2, median)},
                                 error = function(cond){
                                     print('fuu')
                                     return(group2)
                                 })

        g19 = groupAverage1 < 9.5 & groupAverage1 > 8
        g16 = groupAverage1  < 6
        g29 = groupAverage2 < 9.5 & groupAverage2 > 8
        g26 = groupAverage2 < 6



        tempGroupAv2 = vector(length = length(groupAverage2))

        tempGroupAv2[g26 & g19] = tryCatch({(apply(group2[, g26 & g19], 2, max))},
                                           error = function(cond){
                                               print('I hate you damn it!')
                                               if (is.null(nrow(group2))){
                                                   return(group2[g26 & g19])
                                               }else{ return(max(group2[, g26 & g19]))
                                               }})

        tempGroupAv2[g16 & g29] = tryCatch({(apply(group2[, g16 & g29], 2, min))},
                                           error = function(cond){
                                               print('I hate you damn it!')
                                               if (is.null(nrow(group2))){
                                                   return(group2[g16 & g29])
                                               }else{ return(min(group2[, g26 & g19]))
                                               }})


        #groupAverage1[5124]
        #groupAverage2[5124]


        #groupAverage1[7067]
        #groupAverage2[7067]

        add1 = g19 & g26 & groupAverage1>tempGroupAv2
        add2 = g29 & g16 & tempGroupAv2>groupAverage1


        fold = (groupAverage1 - groupAverage2)
        chosen =  which({(fold >= (log(f)/log(2))) & !(g19 & g26) } | {(fold <= log(1/f)/log(2)) &  !(g29 & g16)}| add1 | add2)
        return(
            data.frame(index = chosen, foldChange = fold[chosen])
        )
    }

    #gene index is the index of the gene in exprData, groupInfos are
    giveSilhouette = function(daGeneIndex, groupInfo1, groupInfo2){
        clustering = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo1)+1)
        clustering = clustering[1:nrow(design) %in% c(groupInfo1, groupInfo2)]
        data = (exprData[ (1:nrow(design) %in% c(groupInfo1, groupInfo2)),  daGeneIndex])
        cluster = list(clustering = clustering, data = data)
        silo = silhouette(cluster,dist(data))
        return(mean(silo[,3]))
    }

    #############
    design = read.table(designLoc,header=T,sep='\t')

    allDataPre = read.csv(exprLoc, header = T)
    geneData = allDataPre[,1:3]
    exprData = allDataPre[,4:ncol(allDataPre)]
    exprData = exprData
    design = design[match(colnames(exprData),make.names(design$sampleName),),]

    exprData = t(exprData)

    # deal with region stuff ----
    regions =
        trimNAs(
            trimElement(
                unique(
                    unlist(
                        strsplit(as.character(design[,regionNames]),',')))
                ,c('ALL','All','all')))
    regionBased = expand.grid(groupNames, regions)
    regionGroups = vector(mode = 'list', length = nrow(regionBased))
    names(regionGroups) = paste0(regionBased$Var2,'_',regionBased$Var1)

    for (i in 1:nrow(regionBased)){
        regionGroups[[i]] = design[,as.character(regionBased$Var1[i])]
        regionGroups[[i]][!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),design[,regionNames])] = NA
    }

    design = cbind(design,regionGroups)
    groupNames = c(groupNames, names(regionGroups))

    # get replicate means -----
    # a terrible way to preallocate
    newExpr = exprData[1:length(unique(design$originalIndex)),]
    indexes = unique(design$originalIndex)
    for (i in 1:length(indexes)){
        newExpr[i, ] = apply(exprData[design$originalIndex == indexes[i],], 2,mean)
    }

    newDesign = design[match(indexes,design$originalIndex),]




    nameGroups = vector(mode = 'list', length = len(groupNames))


    names(nameGroups) = c(groupNames)

    for (i in 1:len(groupNames)){
        nameGroups[[i]] = newDesign[,groupNames[i]]
    }



    nameGroups = nameGroups[unlist(lapply(lapply(lapply(nameGroups,unique),trimNAs),length)) > 1]

    stepi = 1
    for (i in nameGroups){
        groupNames = trimNAs(unique(i))
        realGroups = vector(mode = 'list', length = length(groupNames))
        names(realGroups) = groupNames
        for (j in 1:length(groupNames)){
            realGroups[[j]] = which(i == groupNames[j])
        }
        groupAverages = list()

        #take average of every group, tryCatch is for groups with a single member
        for (j in realGroups){
            groupAverage = tryCatch({apply(newExpr[j,], 2, mean)}, #if by itself just output it. but it shouldnt be by itself
                                    error = function(cond){
                                        return(newExpr[j,])
                                        print('something\'s weird, a single replicate for a whole group?')
                                    })
            groupAverages = c(groupAverages, list(groupAverage))
        }

        names(groupAverages)= groupNames
        groupAverages = t(as.data.frame(groupAverages))

        dir.create(outLoc, showWarnings = F)
        dir.create(paste0(outLoc,'/Marker'), showWarnings = F)
        dir.create(paste0(outLoc  , '/', names(nameGroups)[stepi] , '/'), showWarnings = F)
        dir.create(paste0(outLoc , '/Marker/' , names(nameGroups)[stepi] , '/'), showWarnings = F)

        for (j in 1:nrow(groupAverages)){
            fileName = paste0(outLoc  , '/', names(nameGroups)[stepi], '/',  names(realGroups)[j])
            fileName2 = paste0(outLoc , '/Marker/' , names(nameGroups)[stepi] , '/' , names(realGroups)[j])

            #find markers
            isMarker = vector(length = ncol(groupAverages))
            for (t in 1:ncol(groupAverages)){
                isMarker[t] = all(groupAverages[-j, t] + log(10, base=2) < groupAverages[j,t])
            }
            fMarker = data.frame(geneData$Gene.Symbol[isMarker], groupAverages[j,isMarker], tryCatch({apply(groupAverages[-j,isMarker],2,max)}, error = function(e){max(groupAverages[-j,isMarker])}))
            fChange = foldChange(groupAverages[j, ], groupAverages[-j,] )
            fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
            fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]

            #silhouette
            groupInfo1 = which(design[,names(nameGroups)[stepi]] == names(realGroups)[j])
            groupInfo2 = which(design[,names(nameGroups)[stepi]] != names(realGroups)[j] & !is.na(design[,names(nameGroups)[stepi]]))

            silo = vector(length = nrow(fChangePrint))
            for (t in 1:nrow(fChangePrint)){
                silo[t] = giveSilhouette(which(geneData$Gene.Symbol == fChangePrint$geneNames[t]),
                                         groupInfo1,
                                         groupInfo2)
            }

            fChangePrint = cbind(fChangePrint, silo)

            print(fileName)
            # print(i)
            write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
            write.table(fMarker, quote = F, row.names = F, col.names = F, fileName2)

        }
        stepi = stepi + 1
    }

}

