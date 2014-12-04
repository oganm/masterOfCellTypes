GSMassign = function(fileIn,fileOut){
    eval( expr = parse( text = getURL(
        "https://raw.githubusercontent.com/oganm/toSource/master/GEO.r",
        ssl.verifypeer=FALSE) ))

    des = read.table(fileIn,header=T,sep='\t')
    gsms = mapply(gsmFind,des$GSE,des$gsmRegex)

    des$GSM = sapply(gsms,paste,collapse=',')
    write.table(des, fileOut, row.names=FALSE, sep = '\t', quote=F)
}