require(RCurl)
# require(igraph)
require(foreach)
require(doMC)
require(parallel)
require(corpcor)
cores = 4
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)

# avoid parallelization when storing coexpression networks.

eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

# temporary. till migration to run.R ------
coexpLoc = 'Data/humanBootstrap'
genesOut = 'Data/RotSel/Relax'
coexSampOut = 'Data/humanBootGenes'
markerChainOut = 'Data/humanMainNetwork'
filter=T
#####
#markerChain = function(coexpLoc, genesOut, coexSampOut, filter = T){

allGenLocs = list.dirs(genesOut)
allGenLocs = allGenLocs[-1]
geneLists = lapply(allGenLocs, puristOut)
names(geneLists) = basename(allGenLocs)

humanExprs = list.files(coexpLoc, full.names = T)
humanExpr = read.csv(humanExprs[1],
                     header = T ,
                     stringsAsFactors = F)
# this will be used to trim humanExprs and itself converted to humanGene
humanGeneF = humanExpr[,1:3]
source('homologene.R')
fullOrtho = human2mouse(humanGeneF$Gene_Symbol)
source('puristOut.R')

for (i in 3:len(humanExprs)){
    print(i)
    
    humanExpr=read.csv(humanExprs[i],
                       header = T ,
                       stringsAsFactors = F)[,4:ncol(humanExpr)]
    humanExpr = humanExpr[!is.na(humanGeneF$Gene_Symbol),]
    humanGene = humanGeneF[!is.na(humanGeneF$Gene_Symbol),]
    print('file read')
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
        print('filtering Complete')
    }
    #select most variable representative of a gene
    
    dupGenes = names(table(humanGene$Gene_Symbol))[table(humanGene$Gene_Symbol)>1]
    
    dupExpr = humanExpr[humanGene$Gene_Symbol %in% dupGenes,]
    dupGenes = humanGene[humanGene$Gene_Symbol %in% dupGenes,]
    
    newExprData <<- foreach (i = 1:length(unique(dupGenes$Gene_Symbol)), .combine=rbind) %dopar% {
        indexes = which(dupGenes$Gene_Symbol %in% unique(dupGenes$Gene_Symbol)[i])
        groupData = dupExpr[indexes,]
        chosen = which.max(apply(groupData,1,var))
        as.double(groupData[chosen,])
    }
    newGeneData = dupGenes[match(unique(dupGenes$Gene_Symbol),dupGenes$Gene_Symbol),]
    colnames(newExprData) = names(humanExpr)
    humanExpr =  humanExpr[!humanGene$Gene_Symbol %in% dupGenes$Gene_Symbol,]
    humanGene = humanGene[!humanGene$Gene_Symbol %in% dupGenes$Gene_Symbol,]
    humanExpr = rbind(humanExpr,newExprData)
    humanGene = rbind(humanGene,newGeneData)
    print('duplicated genes removed')
    if (i>1){
        rm(cors)
        rm(ranks)
        rm(humanCor)
        rm(relCor)
        rm(relRank)
        print('remnants deleted')
    }
    humanCor = cor(t(humanExpr))
    colnames(humanCor) = humanGene$Gene_Symbol
    rownames(humanCor) = humanGene$Gene_Symbol
    cors = sm2vec(humanCor)
    ranks = rank(-cors)
    ranks = vec2sm(ranks)
    colnames(ranks) = humanGene$Gene_Symbol
    rownames(ranks) = humanGene$Gene_Symbol
    tresh = quantile(cors[cors>0],probs = c(0.99))
    
    print('ranking and shit is done')
    
    
    for (j in 1:len(geneLists)){
        
        dir.create(paste0(markerChainOut,'/mainNetwork/',i,'/',names(geneLists)[j]),
                   showWarnings=F,recursive=T)
        dir.create(paste0(markerChainOut,'/ranks/',i,'/',names(geneLists)[j]), 
                   showWarnings=F,recursive=T)
        dir.create(paste0(markerChainOut,'/subsetCor/',i,'/',names(geneLists)[j]),
                   showWarnings=F,recursive=T)
        
        for (k in 1:len(geneLists[[j]])){
            
            common =  fullOrtho$humanGene[fullOrtho$mouseGene %in% geneLists[[j]][[k]]]
            
 
            subsetCor = humanCor[,humanGene$Gene_Symbol %in% common]
            
            if (is.null(dim(subsetCor))){
                subsetCor = as.data.frame(subsetCor)
                colnames(subsetCor) = common
            }
            
            
            write.csv((humanCor[,humanGene$Gene_Symbol %in% common]>tresh),
                      file = paste0(markerChainOut,'/subsetCor/',i,'/',
                                    names(geneLists)[j],'/',
                                    names(geneLists[[j]])[k]),
                      col.names=T,row.names=T)
            
            write.csv(humanCor[humanGene$Gene_Symbol %in% common,humanGene$Gene_Symbol %in% common]>tresh,
                      file = paste0(markerChainOut,'/mainNetwork/',i,'/',
                                    names(geneLists)[j],'/',
                                    names(geneLists[[j]])[k]),
                      col.names=T,row.names=T)
            
            write.csv(ranks[,humanGene$Gene_Symbol %in% common],
                      file = paste0(markerChainOut,'/ranks/',i,'/',
                                    names(geneLists)[j],'/',
                                    names(geneLists[[j]])[k]),
                      row.names=T,col.names=T)
            
        }    
        
    }
}
#}
