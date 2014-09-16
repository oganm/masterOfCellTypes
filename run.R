
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

# file names ----
finalExp = 'finalExp'
qnormExp= 'qnormExp'

#  for heatMap ----
heatFile = 'images/heatmap.png'
heatProps = c('ourNaming',
              'Reference',
              'Age.of.Mouse..postnatal.day.')
heatGenes = 'ourNaming'
heatPalette = colorRampPalette(c("darkred",'white', "blue"))(n = 1000)
# long coloring var for heatmap. modify to your heart's content ----
# extra cell types will not cause problems if in list
namingColors = c(method = 'direct',
                 legendLoc = 'bottomleft',
                 cex = 0.6,
                 Oligo = 'darkgreen',
                 Bergmann = 'palegreen',
                 MotorCholin = 'red',
                 Cholin = 'darkorange',
                 Spiny = 'blanchedalmond',
                 Gluta = 'slategray',
                 Basket = 'mediumpurple4',
                 Golgi = 'orchid',
                 Pyramidal = 'turquoise',
                 Purkinje = 'purple',
                 Inter = 'pink',
                 Granule = 'thistle',
                 Microglia = 'white',
                 Gaba = 'firebrick4',
                 Astroctyte = 'yellow',
                 GabaPV = 'firebrick2',
                 Stem = 'blue' ,
                 Ependymal = 'orange',
                 Serotonergic = 'darkolivegreen',
                 Hyprocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0')


ageColors = c(method = 'gradient', lo = 'blue', hi = 'red', gradFine = 10000,
              legendLoc = 'bottomright',
              cex = 0.6)

refColors = c(method = 'def', legendLoc = 'none')
# names are not important
heatColors = list(ourNaming = namingColors, Reference = refColors, Age.of.Mouse..postnatal.day. = ageColors)

# dependencies ------
# install.packages('reshape')

# normalization ----
if (skipNorm == F){
source('readDesignMergeCel.R')
readDesignMergeCel(desFile, namingCol, celRegex, celDir,tinyChip, outFolder)

source('quantileNormalize.R')
quantileNorm(paste0(outFolder,'/rmaExp'),
             paste0(outFolder,'/',qnormExp))

source('mostVariable.R')
mostVariable(paste0(outFolder,'/',qnormExp),
             paste0(outFolder,'/',finalExp))
}

if (skipNorm == T){
    source('readDesignMergeCel.R')
    meltDesign(desFile, namingCol, celRegex, paste0(outFolder,'/',finalExp), paste0(outFolder,'/meltedDesign'))
}


source('sexFind.R')
sexFind(paste0(outFolder,'/meltedDesign'),
        paste0(outFolder,'/meltedDesign'),
        paste0(outFolder,'/',finalExp))

# SHREEJOY COMMENT OUT AFTER THIS

# gene selection -----
source('geneSelect.R')
geneSelect(paste0(outFolder,'/meltedDesign'),
           paste0(outFolder,'/',finalExp),
           geneOut,
           groupNames,
           regionNames)


source('microglialException.R')
# intersecting with genes from
# http://www.nature.com/neuro/journal/v17/n1/pdf/nn.3599.pdf
microglialException(geneOut)


# heatmap ----
# might migrate into it's own file if more complicated stuff is required
genes = vector()
for (i in heatGenes){
    filenames = list.files(paste0(geneOut,'/',i),include.dirs = FALSE)
    fileContents = lapply(paste0(geneOut,'/',i,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    for (i in 1:length(fileContents)){
        geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0])
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    genes = c(genes, unlist(puristList))
}


source('heatUp.R')
heatUp(paste0(outFolder,'/',finalExp),
       paste0(outFolder,'/meltedDesign'),
       geneOut,
       heatFile,
       heatProps,
       heatColors,
       genes,
       heatPalette)


heatUp = function(expLoc, designLoc, geneOut, heatFile, heatProps, heatColors, geneList){


# calculates specificity index as describe in
# https://www.landesbioscience.com/journals/systemsbiomedicine/article/25630/
# takes too much time!
# source('specIndex.R')
# specIndex(paste0(outFolder,'/meltedDesign'),
#          paste0(outFolder,'/finalExp'),
#          geneOutIndex,
#          groupNames)


