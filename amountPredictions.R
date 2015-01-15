source('puristOut.R')
puristList = puristOut('Data/RotSel/Relax/Cortex_CellType/')

# for human regions ------

# filtering shit for human data
softFile = read.design('Data/hugeHumanSoft.tsv')
softExpr = read.exp('Data/humanBootstrap/1')
medExpr = median(unlist(softExpr[4:len(softExpr)]))
medVar = median(apply(softExpr[4:len(softExpr)],1,var))
keep = apply(humanExpr,1,function(row){
    mean(row)>medExpr | var(row)>medVar
})
softExpr = softExpr[keep,]

softExpr = mostVariable(softExpr,'Gene_Symbol')


names(softExpr)[4:len(softExpr)] = sub('_.*','',names(softExpr)[4:len(softExpr)] )

groups = softFile$Region[match( names(softExpr)[4:len(softExpr)], softFile$GSM)]

estimates = cellTypeEstimate(softExpr,genes = puristList, geneColName='Gene_Symbol',
                             tableOut=paste0("Data/Estimates/Cortex-White/",names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/Cortex-White/",'indivExp ',names(puristList),'.svg'),
                             groups = groups)

estimates = lapply(estimates,function(x){x[groups %in% c('frontal cortex', 'white matter')]})

groups = groups[groups %in%c('frontal cortex', 'white matter')]

plotEstimates(estimates,groups,
              paste0("Data/Estimates/Cortex-White/",names(puristList),'.svg'))

# for bipolar and scz----
bipolExp = read.exp('Data/Bipol.csv')
bipolDes = read.design('Data/BipolDes.tsv')
puristList = puristOut('Data/RotSel/Relax/Cortex_CellType/')
estimates = cellTypeEstimate(bipolExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipol/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipol/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolDes$disease_state2#,
                             #controlBased='Cont'
                             )
plotEstimates(estimates,bipolDes$disease_state2,
              paste0("Data/Estimates/humanBipol/",names(puristList),'.svg'))

bipolSczExp =  read.exp('Data/BipolScz.csv')
bipolSczDes = read.design('Data/BipolSczDes.tsv')
estimates = cellTypeEstimate(bipolSczExp,genes= puristList,
                             geneColName="Gene.Symbol",
                             tableOut = paste0('Data/Estimates/humanBipolScz/', names(puristList),'.tsv'),
                             indivGenePlot = paste0("Data/Estimates/humanBipolScz/",'indivExp ',names(puristList),'.svg'),
                             groups=bipolSczDes$disease_state,
                             controlBased = 'Cont')
plotEstimates(estimates,bipolSczDes$disease_state,
              paste0("Data/Estimates/humanBipolScz/",names(puristList),'.svg'))
