rotateCheck = function(rotationOut){
    require(RCurl)
    eval( expr = parse( text = getURL(
        "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
        ssl.verifypeer=FALSE) ))

    dirFols = list.dirs(rotationOut, recursive = F)
    if (any(grepl('Confidence',dirFols))){
        dirFols = dirFols[-which(grepl('Confidence',dirFols))]
    }

    loopAround = list.dirs(dirFols[1],full.names = F)
    loopAround = loopAround [-which(loopAround %in% c('Relax','Marker',''))]
    dir.create(paste0(rotationOut,'/Confidence'), showWarnings = F)

    for (i in loopAround){
        dir.create(paste0(rotationOut,'/Confidence/',i),recursive = T, showWarnings = F)
        files = list.files(paste0(dirFols[1],'/',i))
        fileOuts = vector(mode = 'list', length = len(files))
        for (j in files){
            genes = vector()
            for (k in dirFols){
                daFile = tryCatch({read.table(paste0(k,'/',i,'/',j) ,header=F)},
                                  error=function(e){
                                    daFile=data.frame()
                                  })
                genes = c(genes, as.character(daFile$V1))
            }
            geneCounts = table(genes)
            confidence = geneCounts/len(dirFols)

            write.table(as.df(confidence), file = paste0(rotationOut,'/Confidence/',i,'/',j), row.names = F, quote=F, sep='\t')
        }
    }


}