setUpHomologene = function(){
    download.file(url = "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", destfile = 'Data/homologene.data')
    homologene = read.table('Data/homologene.data',sep ='\t',quote='')
    names(homologene) = c('HID','Taxonomy','Gene.ID','Gene.Symbol','Protein.GI','Protein.Accession')
    homoHuman = subset(homologene, Taxonomy == 9606)
    homoHuman = homoHuman[,c('HID','Gene.Symbol')]
    homoMouse = subset(homologene, Taxonomy == 10090)
    rm(homologene)
    homoMouse = homoMouse[,c('HID','Gene.Symbol')]
    names(homoHuman)[2]='humanGene'
    names(homoMouse)[2]='mouseGene'
    homoFile = merge(homoHuman,homoMouse)
    rm('homoHuman','homoMouse')
    system('rm Data/homologene.data')
    write.table(homoFile,file='Data/homologene.tsv', quote=F,sep = '\t',col.names=T,row.names=F)
}


mouse2human = function(genes){
    homolo = read.table('Data/homologene.tsv',header=T,sep='\t')
    homolo = homolo[homolo$mouseGene %in% genes,]
    return(homolo[,c('mouseGene','humanGene')])
}
human2mouse = function(genes){
    homolo = read.table('Data/homologene.tsv',header=T,sep='\t')
    homolo = homolo[homolo$humanGene %in% genes,]
    return(homolo[,c('humanGene','mouseGene')])
}