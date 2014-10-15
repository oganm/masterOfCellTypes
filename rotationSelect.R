rotFold = 'Data/Rotation/Confidence'
rotSelOut = 'Data/RotSel'

require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

tryRead = function (x){
    out = tryCatch({
        read.table(x, stringsAsFactors=F)
    }, error = function(e){
        data.frame()
    })
}

rotationSelect = function(rotFold,geneOut, rotSelOut){
    # in this version of R full.names = F does not work. Shouldnt cause problems in the newer one
    allDirs = gsub(list.dirs(rotFold, full.names=F)[1],'',list.dirs(rotFold, full.names=F))
    allDirs = allDirs[-which(allDirs %in% c('/Marker','/Relax',''))]
    for (i in 1:len(allDirs)){
        daFiles = list.files(paste0(geneOut, allDirs[i]))
        rotationResults = lapply(paste0(rotFold,allDirs[i],'/',daFiles),tryRead)
        names(rotationResults) = daFiles
        foldResults = lapply(paste0(geneOut,allDirs[i],'/',daFiles),tryRead)
        names(foldResults) = daFiles
        
        dir.create(paste0(rotSelOut, allDirs[i]),
                   recursive = T,showWarnings = F)
        
        for (j in 1:len(foldResults)){
            if (all(dim(foldResults[[j]])==c(0,0))){
                if (!all(dim(rotationResults[[j]])==c(0,0))){
                    if (len(rotationResults[[j]][rotationResults[[j]][,2]>0.95,1])>0){
                        print('hmm very interesting...')
                        print(i)
                        print(j)
                    }
                }
                file.create(paste0(rotSelOut, allDirs[i], '/', daFiles[j]))
            } else {
                foldResults[[j]][,1][!foldResults[[j]][,1] %in% rotationResults[[j]][rotationResults[[j]][,2]>0.95,1]]
                
                interesting = rotationResults[[j]][,1] [!rotationResults[[j]][rotationResults[[j]][,2]>0.95,1] %in% foldResults[[j]][,1]]
                if (len(interesting>0)){
                    print('hmm interesting...')
                    print(i)
                    print(j)
                }
                
                inter = intersect(foldResults[[j]][,1], rotationResults[[j]][rotationResults[[j]][,2]>0.95,1])
                
                write.table(foldResults[[j]][foldResults[[j]][,1] %in% inter,],
                            quote = F, row.names = F, col.names = F, 
                            file =paste0(rotSelOut, allDirs[i], '/', daFiles[j]) )
                
            }
        }
    }    
}


# rotationOutput = puristOut('Data/RotSel/Relax/CellType')
# normalOutput = puristOut('Data/Fold/Relax/CellType')
# deLaGenes = list(normal = unlist(normalOutput), rotation = unlist(rotationOutput))
# a = venn.diagram(deLaGenes, filename= NULL, fill = c('red','blue'))
# grid.draw(a)

