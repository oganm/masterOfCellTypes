# for loading and normalization -----
desFile='Data/Design.xls'
xls = TRUE
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'MainName'
tinyChip = 'mouse430a2.db'
skipNorm = T
# for gene selection -----
geneOut = 'Data/Fold'
geneOutIndex = 'Data/Index'
groupNames = c('CellType' , 'FullCellType', 'GabaDeep', 'forContanim')
contanimName = 'forContanim'
regionNames = 'Region'
rotationOut = 'Data/Rotation'

# contamination indeces ----
defMarkers = 'Data/defMarkers.csv'
apContOut = 'Data/UsedMarkers'
# file names ----
finalExp = 'finalExp'
qnormExp= 'qnormExp'

#  for heatMap ----
heatFile = 'images/heatmap.png'
heatProps = c('CellType',
              'MajorType',
              'Reference',
              'Age')
heatGenes = 'CellType'
heatPalette = colorRampPalette(c("darkred",'white', "blue"))(n = 1000)
# long coloring var for heatmap. modify to your heart's content ----
# extra cell types will not cause problems if in list
namingColors = c(method = 'direct',
                 legendLoc = 'bottomleft',
                 cex = 1,
                 Oligo = 'darkgreen',
                 Bergmann = 'palegreen',
                 MotorCholin = 'darkorange4',
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
                 Astrocyte = 'yellow',
                 GabaPV = 'firebrick2',
                 Stem = 'blue' ,
                 Ependymal = 'orange',
                 Serotonergic = 'darkolivegreen',
                 Hyprocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0')


ageColors = c(method = 'gradFactor', lo = 'blue', hi = 'red',
              legendLoc = 'bottomright',factorOrd = NA,
              cex = 1)

typeColors = c(method = 'def', legendLoc = 'none')

refColors = c(method = 'def', legendLoc = 'topleft', cex = 1)
# names are not important
heatColors = list(CellType = namingColors, typeColors = typeColors, Reference = refColors, Age = ageColors)

# dependencies (not complete. will migrate and see) ------
# install.packages('reshape')

# convert to tab delimeted ----
if (xls == TRUE){
    system(paste0('libreoffice --headless --convert-to csv ',desFile))

    splPath = strsplit(desFile,'/')[[1]]
    splPath = splPath[len(splPath)]

    write.table(read.csv(paste0(substr(splPath,1,nchar(splPath)-4),'.csv')), sep = '\t', quote = F, row.names = F,
                file = paste0(substr(desFile,1,nchar(desFile)-4),'.tdf'))
    file.remove(paste0(substr(splPath,1,nchar(splPath)-4),'.csv'))
    desFile = paste0(substr(desFile,1,nchar(desFile)-4),'.tdf')
}


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
microglialException(paste0(geneOut,'/Relax'))
microglialException(paste0(geneOut,'/Marker'))

system('beep')

# contamination -----
source('appendCont.R')
appendCont(defMarkers,paste0(outFolder,'/meltedDesign'),
           paste0(geneOut, '/Marker'), contanimName, apContOut)

source('contaminate.R')
contamination(paste0(outFolder,'/meltedDesign'),
           paste0(outFolder,'/',finalExp),
           apContOut,
           paste0(outFolder,'/meltedDesign'))


source('sampRotate.R')
for (i in 1:10){
    sampRotate(paste0(outFolder,'/meltedDesign'),
               paste0(outFolder,'/',finalExp),
               paste0(rotationOut,'/',i),
               groupNames,
               regionNames)
}


# heatmap ----
# might migrate into it's own file if more complicated stuff is required
expandedHeat = vector()
for (i in heatGenes){
    expandedHeat = c(expandedHeat,
                     list.files(paste0(geneOut,'/Relax/'),include.dirs = T) [grepl(i,list.files(paste0(geneOut,'/Relax/'),include.dirs = T))]
                     )
}
expandedHeat = expandedHeat[-5]

genes = vector()
for (i in expandedHeat){
    filenames = list.files(paste0(geneOut,'/Relax/',i),include.dirs = FALSE)
    fileContents = lapply(paste0(geneOut,'/Relax/',i,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    for (j in 1:length(fileContents)){
        geneList[[j]] = as.character(fileContents[[j]]$V1[(as.numeric(as.character(fileContents[[j]]$V3))>0.5)&(as.numeric(as.character(fileContents[[j]]$V2))>0)])
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    genes = c(genes, unlist(puristList))
}
genes= unique(genes)

source('heatUp.R')
heatUp(paste0(outFolder,'/',finalExp),
       paste0(outFolder,'/meltedDesign'),
       geneOut,
       heatFile,
       heatProps,
       heatColors,
       heatPalette,
       genes,
       2000,
       2000
       )

genes = vector()
for (i in expandedHeat){
    filenames = list.files(paste0(geneOut,'/Relax/',i),include.dirs = FALSE)
    fileContents = lapply(paste0(geneOut,'/Relax/',i,'/', filenames), read.table)
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = filenames
    for (j in 1:length(fileContents)){
        geneList[[j]] = as.character(fileContents[[j]]$V1[as.numeric(as.character(fileContents[[j]]$V3))*as.numeric(as.character(fileContents[[j]]$V2))>2])
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    genes = c(genes, unlist(puristList))
}

genes = unique(genes)

heatUp(paste0(outFolder,'/',finalExp),
       paste0(outFolder,'/meltedDesign'),
       geneOut,
       'images/heatmapMult',
       heatProps,
       heatColors,
       heatPalette,
       genes,
       2000,
       2000
)


heatUp(paste0(outFolder,'/',finalExp),
       paste0(outFolder,'/meltedDesign'),
       geneOut,
       'images/noFilter.png',
       heatProps,
       heatColors,
       heatPalette,
       NA
       )



# calculates specificity index as describe in
# https://www.landesbioscience.com/journals/systemsbiomedicine/article/25630/
# takes too much time!
# source('specIndex.R')
# specIndex(paste0(outFolder,'/meltedDesign'),
#          paste0(outFolder,'/finalExp'),
#          geneOutIndex,
#          groupNames)

system('while true; do beep; sleep 2; done')
