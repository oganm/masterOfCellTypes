require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
require(foreach)
require(doMC)
require(parallel)
require(corpcor)
cores = 4
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)


finalChainOut="Data/CoexpGenes"
finalChain = function(markerChainOut,finalChainOut){
    
    sumNetwork = function(network1,network2){
        rows1 = rownames(network1)
        rows2 = rownames(network2)
        newRows = rows2[!rows2 %in% rows1]
        cols1 = colnames(network1)
        cols2 = colnames(network2)
        newCols = cols2[!cols2 %in% cols1]
        newNetwork = matrix(rep(0,(len(rows1)+len(newRows))*
                                        (len(cols1)+len(newCols))),
                                ncol = len(cols1)+len(newCols))
        rownames(newNetwork) = c(rows1,newRows)
        colnames(newNetwork) = c(cols1,newCols)
        newNetwork[rows1,cols1]=network1
        newNetwork[rows2,cols2] = newNetwork[rows2,cols2]+as.matrix(network2)
        
        return(newNetwork)
    }
    
    toLoop = list.dirs(paste0(markerChainOut,'/subsetCor/1'))[-1]
    toLoop = basename(toLoop)
    allFiles = list.dirs(paste0(markerChainOut,'/subsetCor'))
    
    foreach (i = toLoop) %dopar%{
        toComb = allFiles[grepl(i,allFiles)]
        cellTypeLoop = list.files(toComb[2])
        finalNetwork = matrix()
        for (j in cellTypeLoop){
            for (k in toComb){
                toAdd = read.csv(paste0(k,'/',j),header=T, row.names=1)
                finalNetwork = sumNetwork(finalNetwork,toAdd)
            }
        }
        
    }
    
    
    
}