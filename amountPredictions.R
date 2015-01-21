require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

# for human regions ------
puristList = puristOut('Data/RotSel/Relax/Cortex_GabaDeep/')
# filtering shit for human data
softFile = read.design('Data/hugeHumanSoft.tsv')
softExpr = read.exp('Data/HumanRegionExpr/cortex-white')

names(softExpr)[4:len(softExpr)] = sub('_.*','',names(softExpr)[4:len(softExpr)] )
# determine groups based on sample names from soft file
groups = softFile$Region[match( names(softExpr)[4:len(softExpr)], softFile$GSM)]


medExpr = median(unlist(softExpr[4:len(softExpr)]))
keep = apply(softExpr[4:len(softExpr)],1,function(row){
    return(mean(row[groups %in% unique(groups)[1]])>medExpr & mean(row[groups %in% unique(groups)[2]])>medExpr)
})

mostVarSoft = mostVariable(softExpr,'Gene_Symbol')




fullEstimate(mostVarSoft,
             genes=puristList,
             geneColName="Gene_Symbol",
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

fullEstimate(bipolExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolDes$disease_state2,
             outDir='Data/Estimates/humanBipol/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased='Cont')


bipolSczExp =  read.exp('Data/BipolData/BipolScz.csv')
bipolSczDes = read.design('Data/BipolData/BipolSczDes.tsv')

fullEstimate(bipolSczExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolSczDes$disease_state,
             outDir='Data/Estimates/humanBipolScz/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased='Cont')


# scz bipolar not control based

fullEstimate(bipolExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolDes$disease_state2,
             outDir='Data/Estimates/humanBipolNotControl/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA)

fullEstimate(bipolSczExp,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=bipolSczDes$disease_state,
             outDir='Data/Estimates/humanBipolSczNotControl/',
             seekConsensus=F,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA)
