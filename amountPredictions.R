require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
source('puristOut.R')

# for human regions ------
puristList = puristOut('Data/RotSel/Relax/GabaDeep/')
# filtering shit for human data
softFile = read.design('Data/hugeHumanSoft.tsv')
softExpr = read.exp('Data/HumanRegionExpr/cortex-white')
medExpr = median(unlist(softExpr[4:len(softExpr)]))
medVar = median(apply(softExpr[4:len(softExpr)],1,var))
keep = apply(softExpr[4:len(softExpr)],1,function(row){
    mean(row)>medExpr | var(row)>medVar
})
softExpr = softExpr[keep,]

softExpr = mostVariable(softExpr,'Gene_Symbol')


names(softExpr)[4:len(softExpr)] = sub('_.*','',names(softExpr)[4:len(softExpr)] )

groups = softFile$Region[match( names(softExpr)[4:len(softExpr)], softFile$GSM)]



estimates = cellTypeEstimate(softExpr,genes = puristList[names(puristList) %in% c('Astrocyte','Oligo')], geneColName='Gene_Symbol',
                             tableOut=paste0("Data/Estimates/Cortex-White/",
                                             names(puristList[names(puristList) %in% c('Astrocyte','Oligo')]),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/Cortex-White/",
                                                    names(puristList[names(puristList) %in% c('Astrocyte','Oligo')]),' indivExp','.svg'),
                             groups = groups)
estimates$estimates = trimNAs(estimates$estimates)
estimates$groups = trimNAs(estimates$groups)

groupRotations(softExpr, 
               genes=puristList[names(puristList) %in% c('Astrocyte','Oligo')],
               geneColName = 'Gene_Symbol',
               groups=groups,
               outDir='Data/Estimates/Cortex-White/')

plotEstimates(estimates$estimates,estimates$groups,
              paste0("Data/Estimates/Cortex-White/",
                     names(estimates$estimates),'.svg'))

# for bipolar and scz----
bipolExp = read.exp('Data/BipolData/Bipol.csv')
bipolDes = read.design('Data/BipolData/BipolDes.tsv')
puristList = puristOut('Data/RotSel/Relax/Cortex_GabaDeep/')
estimates = cellTypeEstimate(bipolExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipol/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipol/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolDes$disease_state2,
                             controlBased='Cont'
                             )
estimates$estimates = trimNAs(estimates$estimates)
estimates$groups = trimNAs(estimates$groups)
plotEstimates(estimates$estimates,estimates$groups,
              paste0("Data/Estimates/humanBipol/",names(estimates$estimates),'.svg'))

bipolSczExp =  read.exp('Data/BipolData/BipolScz.csv')
bipolSczDes = read.design('Data/BipolData/BipolSczDes.tsv')
estimates = cellTypeEstimate(bipolSczExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipolScz/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipolScz/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolSczDes$disease_state,
                             controlBased = 'Cont')
estimates$estimates = trimNAs(estimates$estimates)
estimates$groups = trimNAs(estimates$groups)
plotEstimates(estimates$estimates,estimates$groups,
              paste0("Data/Estimates/humanBipolScz/",names(estimates$estimates),'.svg'))


# scz bipolar not control based
puristList = puristOut('Data/RotSel/Relax/Cortex_GabaDeep/')
estimates = cellTypeEstimate(bipolExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipolNotControl/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipolNotControl/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolDes$disease_state2
)
estimates$estimates = trimNAs(estimates$estimates)
estimates$groups = trimNAs(estimates$groups)
plotEstimates(estimates$estimates,estimates$groups,
              paste0("Data/Estimates/humanBipol/",names(estimates$estimates),'.svg'))


estimates = cellTypeEstimate(bipolSczExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipolSczNotControl/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipolSczNotControl/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolSczDes$disease_state)
estimates$estimates = trimNAs(estimates$estimates)
estimates$groups = trimNAs(estimates$groups)
plotEstimates(estimates$estimates,estimates$groups,
              paste0("Data/Estimates/humanBipolScz/",names(estimates$estimates),'.svg'))