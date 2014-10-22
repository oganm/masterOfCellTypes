fuckingOrthologues = function(){
    humanAnot = read.csv('Data/Annotations/HG-U133_Plus_2.na34.annot.csv',comment.char='#')
    mouseAnot = read.csv('Data/Annotations/Mouse430A_2.na34.annot.csv', comment.char ='#')
    orthology = read.csv('Data/Annotations/HG-U133_Plus_2.na34.ortholog.csv', comment.char ='#')
    
    orthology = orthology[orthology$Ortholog.Array == 'Mouse430_2',c('Probe.Set.ID','Ortholog.Probe.Set')]
    orthology$Ortholog.Probe.Set = tolower(orthology$Ortholog.Probe.Set)
    humanAnot = humanAnot[,c('Probe.Set.ID', 'Gene.Symbol')]
    names(humanAnot) = c('Probe.Set.ID','humanGene')
    mouseAnot = mouseAnot[ ,c('Probe.Set.ID','Gene.Symbol')]
    names(mouseAnot) = c('Ortholog.Probe.Set','mouseGene')
    finalOrtho=merge(merge(orthology,humanAnot),mouseAnot)
    write.table(unique(finalOrtho[,c('humanGene','mouseGene')]),file = 'Data/orthology.tsv',quote= F
                ,row.names=F,
                sep = '\t')
}