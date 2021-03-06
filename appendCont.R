appendCont = function(defMarkers, designLoc, markerLoc, contanName, outFile){
    # appends new markers from a folder to a list of known markers also add classification to design file
    # consider separating them
    defM = tryCatch({read.table(defMarkers,header=T,sep='\t')},
                    error = function(e){
                        defM=list()
                    })
    defM = lapply(as.list(defM),trimElement,'')
    design = read.table(designLoc,header=T,sep='\t')

    files = list.files(paste0(markerLoc,'/',contanName),include.dirs = FALSE)
    fileContents = lapply(paste0(markerLoc,'/',contanName,'/', files), function(x){
        tryCatch({
            read.table(x)},
            error = function(e){
                NA
            })})

    names(fileContents) = files

    for (i in 1:len(fileContents)){
        defM[[files[i]]] = tryCatch({unique(c(as.character(defM[[files[i]]]),
                                              as.character(fileContents[[i]]$V1)))},
                                    error = function(e){
                                        c(as.character(defM[[files[i]]]),NULL)
                                    })
        if(!len(defM[[files[i]]])==0){
            if (all(names(design) !=  files[[i]])){
                design[files[[i]]] = design[,contanName] == files[[i]]
            }
        } else {
            defM = defM[-which(names(defM)==files[i])]
        }

    }

    write.table(design, file = designLoc, row.names=FALSE,sep = '\t', quote=F)

    maxLen = max(unlist(lapply(defM,len)))
    for (i in 1:len(defM)){
        defM[[i]] = c(as.character(defM[[i]]) ,rep ('', maxLen - len(defM[[i]])))
    }

    write.table(as.data.frame(defM), file = outFile,  row.names=FALSE,sep = '\t', quote=F)
}