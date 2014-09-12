
# for loading and normalization -----
desFile='Design.csv'
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'MainName'
tinyChip = 'mouse430a2.db'
skipNorm = F
# for gene selection -----
geneOut = 'Data/Fold'
geneOutIndex = 'Data/Index'
groupNames = 'ourNaming'
regionNames = 'Region'

# dependencies ------
# install.packages('reshape')

# executing stuff and things -----
if (skipNorm == F){
source('readDesignMergeCel.R')
readDesignMergeCel(desFile, namingCol, celRegex, celDir,tinyChip, outFolder)

source('quantileNormalize.R')
quantileNorm(paste0(outFolder,'/rmaExp'),
             paste0(outFolder,'/qnormExp'))

source('mostVariable.R')
mostVariable(paste0(outFolder,'/qnormExp'),
             paste0(outFolder,'/finalExp'))
}

if (skipNorm == T){
    source('readDesignMergeCel.R')
    meltDesign(desFile, namingCol, celRegex, paste0(outFolder,'/finalExp'), paste0(outFolder,'/meltedDesign'))
}


source('sexFind.R')
sexFind(paste0(outFolder,'/meltedDesign'),
        paste0(outFolder,'/meltedDesign'),
        paste0(outFolder,'/finalExp'))

# SHREEJOY COMMENT OUT AFTER THIS

source('geneSelect.R')
geneSelect(paste0(outFolder,'/meltedDesign'),
           paste0(outFolder,'/finalExp'),
           geneOut,
           groupNames,
           regionNames)


source('microglialException.R')
# intersecting with genes from
# http://www.nature.com/neuro/journal/v17/n1/pdf/nn.3599.pdf
microglialException(groupNames, geneOut)

# calculates specificity index as describe in
# https://www.landesbioscience.com/journals/systemsbiomedicine/article/25630/
# takes too much time!
# source('specIndex.R')
# specIndex(paste0(outFolder,'/meltedDesign'),
#          paste0(outFolder,'/finalExp'),
#          geneOutIndex,
#          groupNames)


