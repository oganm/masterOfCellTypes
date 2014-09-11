
# for loading and normalization -----
desFile='Design.csv'
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'MainName'
tinyChip = 'mouse430a2.db'

# for gene selection -----
geneOut = 'Data/Rest'
groupNames = 'ourNaming'

# dependencies ------
# install.packages('reshape')

# executing stuff and things -----
source('readDesignMergeCel.R')

readDesignMergeCel(desFile, namingCol, celRegex, celDir,tinyChip, outFolder)

source('quantileNormalize.R')
quantileNorm(paste0(outFolder,'/allNormalized'),
             paste0(outFolder,'/quantileNormalized'))

source('mostVariable.R')
mostVariable(paste0(outFolder,'/quantileNormalized'),
             paste0(outFolder,'/mostVariableQuantileNormalized'))


source('sexFind.R')
sexFind(paste0(outFolder,'/normalizedDesign'),
        paste0(outFolder,'/normalizedDesign'),
        paste0(outFolder,'/mostVariableQuantileNormalized'))

source('geneSelect.R')
geneSelect(paste0(outFolder,'/normalizedDesign'),
           paste0(outFolder,'/mostVariableQuantileNormalized'),
           geneOut,
           groupNames)


source('microglialException.R')
microglialException(groupNames, geneOut)




function(designLoc,exprLoc,outLoc,groupNames)
    designLoc= paste0(outFolder,'/normalizedDesign')
exprLoc = paste0(outFolder,'/mostVariableQuantileNormalized')
outLoc = geneOut
groupNames = groupNames

# microglia exception
