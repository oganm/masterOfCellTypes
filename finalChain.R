require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
require(foreach)
require(doMC)
require(parallel)
require(corpcor)
cores = 8
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)


finalChainOut="Data/CoexpNetwork"
markerChainOut = 'Data/humanMainNetwork/'
# finalChain = function(markerChainOut,finalChainOut){
    
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
    allFiles = allFiles[!grepl('76',allFiles)]
    foreach (i = toLoop) %dopar%{
        toComb = allFiles[grepl(i,allFiles)]
        cellTypeLoop = list.files(toComb[1])
        for (j in cellTypeLoop){
            finalNetwork = matrix()
            for (k in toComb){
                tryCatch({
                    toAdd = read.csv(paste0(k,'/',j),header=T, row.names=1)
                    # for some reason rownames are not added to subsetCor files. 
                    # you have to read the cornames from ranking files.
                    # this could be avoided if you weren't such a moron.
                    moronicBaffoon=read.csv(gsub('subsetCor','ranks',paste0(k,'/',j)),header=T, row.names=1)
                    rownames(toAdd) = rownames(moronicBaffoon)
                    finalNetwork = sumNetwork(finalNetwork,toAdd)
                    },
                    error = function(e){
                        print("nothing to see 'ere")
                    })
            }
            finalNetwork = finalNetwork/len(toComb)
            dir.create(paste0(finalChainOut,'/',i),recursive=T,showWarnings=F)
            write.csv(finalNetwork,paste0(finalChainOut,'/',i,'/',j))
            print(i)
        }
        
    }    
#}