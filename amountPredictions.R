require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
source('predictAmount.R')
source('puristOut.R')
source('runVars.R')
require('reshape2')
require(ggplot2)
source('mostVariable.R')
# for human regions ------
puristList = puristOut('Data/RotSel/Relax/Cortex_GabaDeep/')
# filtering shit for human data
softFile = read.design('Data/hugeHumanSoft.tsv')
softExpr = read.exp('Data/HumanRegionExpr/Mixed/cortex-white')

names(softExpr)[4:len(softExpr)] = sub('_.*','',names(softExpr)[4:len(softExpr)] )
# determine groups based on sample names from soft file
list[softGenes,softExp] = sepExpr(softExpr)
groups = softFile$Region[match(names(softExp), softFile$GSM)]


medExpr = median(unlist(softExp))
keep = apply(softExp,1,function(row){
    return(mean(row[groups %in% unique(groups)[1]])>medExpr & mean(row[groups %in% unique(groups)[2]])>medExpr)
})
softExpr=softExpr[keep,]
mostVarSoft = mostVariable(softExpr,'Gene.Symbol')




fullEstimate(mostVarSoft,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=groups,
             outDir='Data/Estimates/Cortex-White/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA)


# for bipolar and scz----
bipolExp = read.exp('Data/BipolData/Bipol.csv')
bipolDes = read.design('Data/BipolData/BipolDes.tsv')
puristList = puristOut('Data/RotSel/Relax/Cortex_GabaDeep/')

bipolExp = bipolExp[,!names(bipolExp) %in% c('BP_18', 'BP_11', 'BP_10', 'Cont_27', 'BP_5')]
bipolDes = bipolDes[!bipolDes$characteristics %in% c('BP_18', 'BP_11', 'BP_10', 'Cont_27', 'BP_5'),]

fullEstimate(bipolExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolDes$disease_state2,
             outDir='Data/Estimates/humanBipol/',
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
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolSczDes$disease_state,
             outDir='Data/Estimates/humanBipolScz/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased='Cont',
             pAdjMethod='none')


# scz bipolar not control based

fullEstimate(bipolExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolDes$disease_state2,
             outDir='Data/Estimates/humanBipolNotControl/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA,
             pAdjMethod='none')

fullEstimate(bipolSczExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolSczDes$disease_state,
             outDir='Data/Estimates/humanBipolSczNotControl/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA,
             pAdjMethod='none')
