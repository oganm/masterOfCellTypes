library(biomaRt)

metaphorDir = 'metaphors/phylomedb.org/metaphors'
metaphorPrep = function(metaphorDir){
    # locate mouse and human files -----
    system(paste0('gunzip ',metaphorDir,'/latest/','species.txt.gz'))
    species = read.table(paste0(metaphorDir,'/latest/','species.txt'),header=T,sep = '\t',quote='"')
    system(paste0('gzip ',metaphorDir,'/latest/','species.txt'))
    
    human = species[species[3] =='Homo sapiens',1]
    mouseCode = species[species[3] =='Mus musculus',2]
    mouse = species[species[3] =='Mus musculus',1]
    humanCode = species[species[3] =='Homo sapiens',2]
    
    
    # matching metaphor IDs to gene symbols -----
    system(paste0('gunzip ',metaphorDir,'/latest/','id_conversion.txt.gz'))
    
    # isolate the mouse and human ones and just take swissport IDs
    idConvert = read.table(paste0(metaphorDir,'/latest/','id_conversion.txt'),header=T,sep = '\t',quote='"',comment.char = '')
    idConvert = idConvert[grepl(paste0('_','((',mouseCode,')|(',humanCode,'))$'),
                                idConvert$protid),]
    idConvert = idConvert[idConvert$db %in% c('swissprot','trembl'),]

    
    system(paste0('gzip ',metaphorDir,'/latest/','id_conversion.txt &'))
    
    # get gene translations
    humanMart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    humanTrans = getBM(attributes = c('uniprot_swissprot','hgnc_symbol'),
                       # just take the human ones. just in case...
                       filters = 'uniprot_swissprot',
                       values = idConvert$X.extid,
                       mart = humanMart)
    humanTrans2 = getBM(attributes = c('uniprot_sptrembl','hgnc_symbol'),
                        # just take the human ones. just in case...
                        filters = 'uniprot_sptrembl',
                        values = idConvert$X.extid,
                        mart = humanMart)
    


    names(humanTrans) = c('X.extid', 'symbol')
    names(humanTrans2) = c('X.extid', 'symbol')
    
    humanTrans=rbind(humanTrans,humanTrans2)
    
    mouseMart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
    
    
    mouseTrans = getBM(attributes = c('uniprot_swissprot','mgi_symbol'),
                       # just take the human ones. just in case...
                       filters = 'uniprot_swissprot', 
                       values = idConvert$X.extid,
                       mart = mouseMart)
    mouseTrans2 = getBM(attributes = c('uniprot_sptrembl','mgi_symbol'),
                       # just take the human ones. just in case...
                       filters = 'uniprot_sptrembl', 
                       values = idConvert$X.extid,
                       mart = mouseMart)
    names(mouseTrans) = c('X.extid','symbol')
    names(mouseTrans2) = c('X.extid','symbol')
    mouseTrans = rbind(mouseTrans,mouseTrans2)
    
    allTrans = rbind(mouseTrans,humanTrans)
    
    mouseTranslate = merge(idConvert, mouseTrans)
    humanTranslate = merge(idConvert, humanTrans)
    symbolSwiss =  merge(idConvert,allTrans)
    
    
    
  
    
    
    # get human orthologues ------
    system(paste0('gunzip ',metaphorDir,'/latest/','orthologs/',human,'.txt.gz'))
    
    humanOrthos =  read.table(paste0(metaphorDir,'/latest/','orthologs/',human,'.txt'),header=T,sep = '\t',quote='"',comment.char = '')
    
    system(paste0('gzip ',metaphorDir,'/latest/','orthologs/',human,'.txt &'))
    
    humanOrthos=humanOrthos[grepl(paste0('_',mouseCode,'$'),humanOrthos$protid2),]
    
    names(humanTranslate)[c(3,4)] = c('protid1','humanGene')

    humanOrthos = merge(humanOrthos, humanTranslate[,3:4])
    names(mouseTranslate)[c(3,4)] = c('protid2', 'mouseGene')
    humanOrthos = merge(humanOrthos, mouseTranslate[,3:4])
    write.table(unique(humanOrthos[,c('humanGene','mouseGene','CS')]), file = 'Data/homology/humanMouse',quote= F
                ,row.names=F,
                sep = '\t')
    
    
    # get mouse orthologues -------
    system(paste0('gunzip ',metaphorDir,'/latest/','orthologs/',mouse,'.txt.gz'))
    mouseOrthos =  read.table(paste0(metaphorDir,'/latest/','orthologs/',mouse,'.txt'),header=T,sep = '\t',quote='"',comment.char = '')
    system(paste0('gzip ',metaphorDir,'/latest/','orthologs/',mouse,'.txt &'))
    mouseOrthos=mouseOrthos[grepl(paste0('_',humanCode,'$'),mouseOrthos$protid2),]
    names(mouseTranslate)[c(3,4)] = c('protid1','mouseGene')
    mouseOrthos = merge(mouseOrthos, mouseTranslate[,3:4])
    names(humanTranslate)[c(3,4)] = c('protid2', 'humanGene')
    mouseOrthos = merge(mouseOrthos, humanTranslate[,3:4])
    write.table(unique(mouseOrthos[,c('mouseGene','humanGene','CS')]), file = 'Data/homology/mouseHuman',quote= F
                ,row.names=F,
                sep = '\t')
    

}

human2mouse = function(geneList){
    table = read.table('Data/homology/humanMouse',header= T, sep = '\t')
    return(unique(table$mouseGene[table$humanGene %in% geneList]))
}

mouse2human = function(geneList){
    table = read.table('Data/homology/mouseHuman',header= T, sep = '\t')
    return(unique(table$humanGene[table$mouseGene %in% geneList]))
}

