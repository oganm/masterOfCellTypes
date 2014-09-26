rotateCheck = function(rotationOut){
    require(RCurl)
    eval( expr = parse( text = getURL(
        "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
        ssl.verifypeer=FALSE) ))

    dirFols = list.dirs(rotationOut, recursive = F)
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
                daFile = read.table(paste0(k,'/',i,'/',j) ,header=F)

            }
        }
    }


}