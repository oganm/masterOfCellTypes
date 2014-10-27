require(RCurl)
require(igraph)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

# temporary. till migration to run.R ------
coexpLoc = 'Data/humanBootstrap'
genesOut = 'Data/Fold/RotSel'
coexSampOut = 'Data/humanBootGenes'
#####
coexpSampleAna = function(coexpLoc, genesOut, coexpHOut, filter = T){
    allGenLocs = list.dirs(genesOut)
    allGenLocs = allGenLocs[-1]
    humanExprs = list.files(coexpLoc, full.names = T)
    humanExpr = read.csv(humanExprs[1],
                         header = T ,
                         stringsAsFactors = F)
    # this will be used to trim humanExprs and itself converted to humanGene
    humanGeneF = humanExpr[,1:3]
    source('homologene.R')
    fullOrtho = human2mouse(humanGeneF$Gene_Symbol)
    source('puristOut.R')
    
    for (i in 1:len(humanExprs)){
        humanExpr=read.csv(humanExprs[i],
                           header = T ,
                           stringsAsFactors = F)[,4:ncol(humanExpr)]
        humanExpr = humanExpr[!is.na(humanGeneF$Gene_Symbol),]
        humanGene = humanGeneF[!is.na(humanGeneF$Gene_Symbol),]
        
        # filtering as described in http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0003911#s3
        # slightly modified to use variance instead of range but the spirit is the same
        if (filter == T){
            medExp = median(unlist(humanExpr))
            medVar = median(apply(humanExpr,1,var))
            keep = apply(humanExpr,1,function(row){
                mean(row)>medExp | var(row)>medVar
            })
            humanExpr = humanExpr[keep,]
            humanGene = humanGene[keep,]
        }
        
        humanCor = cor(t(humanExpr))
        
        cors = humanCor[upper.tri(humanCor)]
        cors = cors[cors>0]
        tresh = quantile(cors,probs = c(0.99))
        
        ranks = rank(-unlist(humanCor),ties.method ='min')
        ranks = matrix(unlist(ranks),ncol=ncol(humanCor))
        
        geneLists = lapply(allGenLocs, puristOut)
        names(geneLists) = basename(allGenLocs)
        
        
        for (j in 1:len(geneLists)){
            dir.create(paste0(coexpHOut,'/list1/',i,'/',
                              names(geneLists[j])),
                       recursive =T, showWarnings = F)
            dir.create(paste0(coexpHOut,'/list0/',i,'/',
                              names(geneLists[j])),
                       recursive =T, showWarnings = F)
            
            dir.create(paste0(coexpHOut,'/ranks/',i,'/',
                              names(geneLists[j])),
                       recursive =T, showWarnings = F)
            dir.create(paste0(coexpHOut,'/used/',i,'/',
                              names(geneLists[j])),
                       recursive =T, showWarnings = F)
            
            for (k in 1:len(geneLists[[j]])){
                common =  fullOrtho$humanGene[fullOrtho$mouseGene %in% geneLists[[j]][[k]]]
                # subset corralation matrix
                relCor = humanCor[,humanGene$Gene_Symbol %in% common]
                
                # adjGraph = graph.adjacency(relCor[humanGene$Gene_Symbol %in% common,]>tresh&relCor[humanGene$Gene_Symbol %in% common,]!=1)
                
                colnames(relCor) = humanGene$Gene_Symbol[humanGene$Gene_Symbol %in% common]
                rownames(relCor) = humanGene$Gene_Symbol
                relRank= ranks[,humanGene$Gene_Symbol %in% common]
                colnames(relRank) = humanGene$Gene_Symbol[humanGene$Gene_Symbol %in% common]
                rownames(relRank) = humanGene$Gene_Symbol
                write.csv(relRank , file =paste0(coexpHOut,'/ranks/',i,'/',
                                                 names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                          quote= F,row.names=T,col.names=T)
                connections = vector(mode = 'list', length = ncol(relCor))
                
                names(connections) = colnames(relCor)
                for (l in 1:ncol(relCor)){
                    connections[[l]] = rownames(relCor)[relCor[,l]>tresh]
                }
                
                
                
                interConnect = sapply(names(connections),findInList,connections)
                
                if (max(sapply(interConnect,len))<=1){
                    file.create(paste0(coexpHOut,'/list1/',i,'/',
                                       names(geneLists[j])))
                    file.create(paste0(coexpHOut,'/list/',i,'/',
                                       names(geneLists[j]),'/',names(geneLists[[j]][k])))
                    file.create(paste0(coexpHOut,'/used/',i,'/',
                                       names(geneLists[j]),'/',names(geneLists[[j]][k])))
                    next 
                }
                
                
                connectionCounts=table(unlist(connections[interConnect[[which.max(sapply(interConnect,len))]]]))
                
                humanList1=names(connectionCounts[connectionCounts>1])
                write.table(humanList1 ,
                            file =paste0(coexpHOut,'/list1/',i,'/',
                                         names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                            quote= F, col.names = F, row.names = F)
                
                humanList = unique(unlist(connections[interConnect[[which.max(sapply(interConnect,len))]]]))
                write.table(humanList ,
                            file =paste0(coexpHOut,'/list/',i,'/',
                                         names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                            quote= F, col.names = F, row.names = F)
                
                
                write.table(names(interConnect)[interConnect[[which.max(sapply(interConnect,len))]]] ,
                            file =paste0(coexpHOut,'/used/',i,'/',
                                         names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                            quote= F, col.names = F, row.names = F)
                
                print(paste(i,j,k))
                # writes the ranks of original marker genes and expanded marker genes.
                # markerRanks=relRank[rownames(relRank) %in% humanList,]
                # write.csv(markerRanks ,
                #           file =paste0(coexpHOut,'/ranks/',i,'/',
                #                        names(geneLists[j]),'/',names(geneLists[[j]][k])) ,
                #           quote= F,row.names=T,col.names=T)
                
                
                
            }
            
            
        }
    }
}
