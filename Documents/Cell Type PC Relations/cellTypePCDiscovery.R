# this aims to check for principle component correlations to metadata and each other for marker genes.
# for microglia, we discover that genes that are differentially expressed in activated microglia 
# mess up the results significantly. which is not nice...
bipolSczExp =  read.exp('Data/BipolData/BipolScz.csv')
bipolSczDes = read.design('Data/BipolData/BipolSczDes.tsv')

bipolSczExp = bipolSczExp[,!names(bipolSczExp) %in% c('BP_23', 'Cont_30', 'BP_12', 'Cont_34')]
bipolSczDes = bipolSczDes[!bipolSczDes$characteristics %in% c('BP_23', 'Cont_30', 'BP_12', 'Cont_34'),]


list[sczGenes,bipolSczExpr] = sepExpr(bipolSczExp)
bipolSczExpr = bipolSczExpr[,bipolSczDes$disease_state %in% 'Cont'] 
bipolSczDes = bipolSczDes[bipolSczDes$disease_state %in% 'Cont',]


bipolExp = read.exp('Data/BipolData/Bipol.csv')
bipolDes = read.design('Data/BipolData/BipolDes.tsv')

bipolExp = bipolExp[,!names(bipolExp) %in% c('BP_18', 'BP_11', 'BP_10', 'Cont_27', 'BP_5')]
bipolDes = bipolDes[!bipolDes$characteristics %in% c('BP_18', 'BP_11', 'BP_10', 'Cont_27', 'BP_5'),]

list[bipolGenes,bipolExpr] = sepExpr(bipolExp)
bipolExpr = bipolExpr[,bipolDes$disease_state %in% 'Cont'] 
bipolDes = bipolDes[bipolDes$disease_state %in% 'Cont',]



commonGenes = intersect(sczGenes$Gene.Symbol,bipolGenes$Gene.Symbol)

bipolSczExpr = bipolSczExpr[match(commonGenes, sczGenes$Gene.Symbol),]
bipolExpr = bipolExpr[match(commonGenes, bipolGenes$Gene.Symbol),]



sczPC = prcomp(t(bipolSczExpr),scale = T )
bpPC = prcomp(t(bipolExpr), scale = T)

toHeat = cbind(sczPC$rotation[,1:4],bpPC$rotation[,1:4])
toHeat= cor(toHeat)
colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
heatmap.2(abs(toHeat),trace='none')

phScz = cor(sczPC$x[,1:4], bipolSczDes$pH)
phBp  = cor(bpPC$x[,1:4], bipolDes$pH)
heatmap.2(abs(cbind(phScz,phBp)),trace='none')

phScz = cor(sczPC$x[,1:4], bipolSczDes$age)
phBp  = cor(bpPC$x[,1:4], bipolDes$age)
heatmap.2(abs(cbind(phScz,phBp)),trace='none')


microSczPC = prcomp(t(bipolSczExpr[commonGenes %in% mouse2human(ourGenes$Microglia)$humanGene,]), scale =T)
microBpPC =  prcomp(t(bipolExpr[commonGenes %in% mouse2human(ourGenes$Microglia)$humanGene,]), scale= T)

toHeat = cbind(microSczPC$rotation[,1:4],microBpPC$rotation[,1:4])
colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
toHeat= cor(toHeat)
heatmap.2(abs(toHeat),trace='none')

dir.create('Documents/Cell Type PC Relations')
for (i in 1:len(ourGenes)){
    tryCatch({
        dir.create(paste0('Documents/Cell Type PC Relations','/',names(ourGenes)[i]))
        microSczPC = prcomp(t(bipolSczExpr[commonGenes %in% mouse2human(ourGenes[[i]])$humanGene,]), scale =T)
        microBpPC =  prcomp(t(bipolExpr[commonGenes %in% mouse2human(ourGenes[[i]])$humanGene,]), scale= T)
        toHeat = cbind(microSczPC$rotation[,1:4],microBpPC$rotation[,1:4])
        colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
        toHeat= cor(toHeat)
        svg(paste0('Documents/Cell Type PC Relations/',names(ourGenes)[i],'/PC cors'))
        heatmap.2(abs(toHeat),trace='none', main = names(ourGenes)[i])
        dev.off()
        phScz = cor(microSczPC$x[,1:4], bipolSczDes$pH)
        phBp  = cor(microBpPC$x[,1:4], bipolDes$pH)
        toHeat = data.frame(SCZ.BP = phScz, BP = phBp)
        svg(paste0('Documents/Cell Type PC Relations/',names(ourGenes)[i],'/pH'))
        heatmap.2(abs(as.matrix(toHeat)),trace='none', cexCol= 1.6 ,main = paste0(names(ourGenes)[i], " pH"), Rowv = F, Colv =F)
        dev.off()
        phScz = cor(microSczPC$x[,1:4], bipolSczDes$age)
        phBp  = cor(microBpPC$x[,1:4], bipolDes$age)
        toHeat = data.frame(SCZ.BP = phScz, BP = phBp)
        svg(paste0('Documents/Cell Type PC Relations/',names(ourGenes)[i],'/age'))
        heatmap.2(abs(as.matrix(toHeat)),trace='none', cexCol= 1.6 ,main = paste0(names(ourGenes)[i], " age"), Rowv = F, Colv =F)
        dev.off()
        phScz = cor(microSczPC$x[,1:4], bipolSczDes$post_mortem)
        phBp  = cor(microBpPC$x[,1:4], bipolDes$post_mortem)
        toHeat = data.frame(SCZ.BP = phScz, BP = phBp)
        svg(paste0('Documents/Cell Type PC Relations/',names(ourGenes)[i],'/postMort'))
        heatmap.2(abs(as.matrix(toHeat)),trace='none', cexCol= 1.6 ,main = paste0(names(ourGenes)[i], " postMort"), Rowv = F, Colv =F)
        dev.off()
        
    } , error =function (e) {
        
    })
}
    
# stimulated genes-------------------
url = "http://bioinf.nl:8080/GOAD/databaseSelectServlet?comparisonID=on&comparisonID=on&fold_change=ALL&p_value=0.01&comparison_ids=22%2C23"

site = getURL(url)

list = str_extract(site, perl('(?<=(var vennData = )).*') )
# separate components
list = unlist(strsplit(unlist(list), '([}]|[\"]),([\"]|[{])'))
# find group names and locations
groupHeads = grep('name',list)
groupNames = str_extract(list[groupHeads],perl('(?<=:[\"]).*'))
# some cleaning
list[grep('data',list)] = str_extract(list[grep('data',list)], perl('(?<=:[[][\"]).*'))
list[grep('}|[]]',list)] =  str_extract(list[grep('}|[]]',list)], perl('^.*?(?=[\"])'))


geneList = vector(mode = 'list', length = length(groupNames))
names(geneList) = groupNames
for (i in 1:length(groupHeads)){
    if (i != length(groupHeads)){
        geneList[[i]] = list[(groupHeads[i]+1):(groupHeads[i+1])-1]
    } else{
        geneList[[i]] =list[(groupHeads[i]+1):length(list)]
    }
}

# check microglia stuff without stimulated genes
toUpper(names(ourGenes)[5]) %in% geneList[[1]]

noFucksGiven  = ourGenes[['Microglia']][!toupper(ourGenes[['Microglia']]) %in% geneList[[1]]]
microSczPC = prcomp(t(bipolSczExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale =T)
microBpPC =  prcomp(t(bipolExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale= T)
toHeat = cbind(microSczPC$rotation[,1:4],microBpPC$rotation[,1:4])
colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
toHeat= cor(toHeat)
heatmap.2(abs(toHeat),trace='none', main = names(ourGenes['Microglia']))


# you re-run the fold gene selection without applyting the microglialException before this
microAll = puristOut('Data/Fold/Relax/Cortex_PyramidalDeep')$Microglia
noFucksGiven = microAll[!toupper(microAll) %in% geneList[[1]]]
microSczPC = prcomp(t(bipolSczExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale =T)
microBpPC =  prcomp(t(bipolExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale= T)
toHeat = cbind(microSczPC$rotation[,1:4],microBpPC$rotation[,1:4])
colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
toHeat= cor(toHeat)

svg('Documents/Cell Type PC Relations/Microglia - ignore stim')
heatmap.2(abs(toHeat),trace='none', main = 'Microglia - ignore stim')
dev.off()

microAll = puristOut('Data/Fold/Relax/Cortex_PyramidalDeep')$Microglia
noFucksGiven = microAll[toupper(microAll) %in% geneList[[1]]]
microSczPC = prcomp(t(bipolSczExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale =T)
microBpPC =  prcomp(t(bipolExpr[commonGenes %in% mouse2human(noFucksGiven)$humanGene,]), scale= T)
toHeat = cbind(microSczPC$rotation[,1:4],microBpPC$rotation[,1:4])
colnames(toHeat) = c('1PC1','1PC2','1PC3','1PC4','2PC1','2PC2','2PC3','2PC4')
toHeat= cor(toHeat)

svg('Documents/Cell Type PC Relations/Microglia - just stim')
heatmap.2(abs(toHeat),trace='none', main = 'Microglia - just stim')
dev.off()

# amount check quick and dirty------
hede = puristOut('Data/Fold/Relax/Cortex_PyramidalDeep')
hede$Microglia = hede$Microglia[!toupper(hede$Microglia) %in% geneList[[1]]]
fullEstimate(bipolExp,
             genes=hede,
             geneColName="Gene.Symbol",
             groups=bipolDes$disease_state2,
             outDir='Data/Estimates/TempBipol/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased='Cont',
             pAdjMethod='none')

bipolSczExp =  read.exp('Data/BipolData/BipolScz.csv')
bipolSczDes = read.design('Data/BipolData/BipolSczDes.tsv')

bipolSczExp = bipolSczExp[,!names(bipolSczExp) %in% c('BP_23', 'Cont_30', 'BP_12', 'Cont_34')]
bipolSczDes = bipolSczDes[!bipolSczDes$characteristics %in% c('BP_23', 'Cont_30', 'BP_12', 'Cont_34'),]

fullEstimate(bipolSczExp,
             genes=hede,
             geneColName="Gene.Symbol",
             groups=bipolSczDes$disease_state,
             outDir='Data/Estimates/TempScz/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased='Cont',
             pAdjMethod='none')


