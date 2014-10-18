require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r",
    ssl.verifypeer=FALSE) ))

# for loading and normalization -----
desFile='Data/Design.tsv'
xls = F
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'Normalize'
tinyChip = 'mouse430a2.db'
skipNorm = F
# for gene selection -----
geneOut = 'Data/Fold'
geneOutIndex = 'Data/Index'
groupNames = c('CellType' , 'FullCellType', 'GabaDeep', 'forContanim')
contanimName = 'forContanim'
regionNames = 'Region'
rotationOut = 'Data/Rotation'
rotSelOut = 'Data/RotSel'
# contamination indeces ----
defMarkers = 'Data/defMarkers.csv'
apContOut = 'Data/UsedMarkers.tsv'

# coexpression stuff
humanDes = 'Data/hugeHumanSoft.tsv'
regionMapping = 'Data/regionMapping.csv'
humanDat = 'Data/HumanExpr'
humanOrtho = 'Data/HT_HG-U133A.na34.ortholog.csv'
coexpHOut = 'Data/HumanMarkerCoexpressed'
# file names ----
finalExp = 'finalExp.csv'
qnormExp= 'qnormExp.csv'



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
                 cex = 1.5,
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
                 CerebGranule = 'thistle',
                 DentateGranule = 'thistle3',
                 Microglia = 'white',
                 Gaba = 'firebrick4',
                 Astrocyte = 'yellow',
                 GabaPV = 'firebrick2',
                 Stem = 'blue' ,
                 Ependymal = 'orange',
                 Serotonergic = 'darkolivegreen',
                 Hyprocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0',
                 Th_positive_LC = 'blueviolet')


ageColors = c(method = 'gradFactor', lo = 'blue', hi = 'red',
              legendLoc = 'bottomright',factorOrd = NA,
              cex = 1.5)

typeColors = c(method = 'def', legendLoc = 'none')

refColors = c(method = 'def', legendLoc = 'topleft', cex = 1.5)
# names are not important
heatColors = list(CellType = namingColors, typeColors = typeColors, Reference = refColors, Age = ageColors)

# for bipol ----
bipolLoc = 'Data/GSE5388_GSE12649_for_Ogan.RData'
bipolOut = 'Data/BipolPC'

# dependencies (not complete. will migrate and see) ------
# install.packages('reshape')

# convert to tab delimeted ----
# does not work in the server
if (xls == TRUE){
    system(paste0('libreoffice --headless --convert-to csv ',desFile))

    splPath = strsplit(desFile,'/')[[1]]
    splPath = splPath[len(splPath)]

    write.table(read.csv(paste0(substr(splPath,1,nchar(splPath)-4),'.csv')), sep = '\t', quote = F, row.names = F,
                file = paste0(substr(desFile,1,nchar(desFile)-4),'.tsv'))
    file.remove(paste0(substr(splPath,1,nchar(splPath)-4),'.csv'))
    desFile = paste0(substr(desFile,1,nchar(desFile)-4),'.tsv')
    #desFile = paste0(substr(desFile,1,nchar(desFile)-4),'.xls')

}


# normalization ----
if (skipNorm == F){
    source('readDesignMergeCel.R')
    readDesignMergeCel(desFile, namingCol, celRegex, celDir,tinyChip, outFolder)


    source('quantileNormalize.R')
    quantileNorm(paste0(outFolder,'/rmaExp.csv'),
                 paste0(outFolder,'/',qnormExp))
    system(paste0('gzip ',outFolder,'/rmaExp.csv'))

    source('mostVariable.R')
    mostVariable(paste0(outFolder,'/',qnormExp),
                 paste0(outFolder,'/',finalExp))
    system(paste0('gzip ',outFolder,'/',qnormExp))
    
}

if (skipNorm == T){
    source('readDesignMergeCel.R')
    meltDesign(desFile, namingCol, celRegex, paste0(outFolder,'/',finalExp), paste0(outFolder,'/meltedDesign.tsv'))
}


source('sexFind.R')
sexFind(paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/',finalExp))




# SHREEJOY COMMENT OUT AFTER THIS

# gene selection -----
source('geneSelect.R')

geneSelect(paste0(outFolder,'/meltedDesign.tsv'),
           paste0(outFolder,'/',finalExp),
           geneOut,
           groupNames,
           regionNames)




source('microglialException.R')
# intersecting with genes from
# http://www.nature.com/neuro/journal/v17/n1/pdf/nn.3599.pdf
microglialException(paste0(geneOut,'/Relax'))
microglialException(paste0(geneOut,'/Marker'))

#system('beep')

# contamination -----
source('appendCont.R')
appendCont(defMarkers,paste0(outFolder,'/meltedDesign.tsv'),
           paste0(geneOut, '/Marker'), contanimName, apContOut)

source('contaminate.R')
contamination(paste0(outFolder,'/meltedDesign.tsv'),
           paste0(outFolder,'/',finalExp),
           apContOut,
           paste0(outFolder,'/meltedDesign.tsv'))


# heatmap -----
source('heatUp.R')
source('heatGeneOut.R')
genes = heatGeneOut(paste0(geneOut,'/Relax/'), heatGenes, 2, T)


heatUp(paste0(outFolder,'/',finalExp),
       paste0(outFolder,'/meltedDesign.tsv'),
       heatFile,
       heatProps,
       heatColors,
       heatPalette,
       genes,
       2000,
       2000)


# sample rotate -----
source('sampRotate.R')

for (i in 1:100){
    sampRotate(paste0(outFolder,'/meltedDesign.tsv'),
               paste0(outFolder,'/',finalExp),
               paste0(rotationOut,'/',i),
               groupNames,
               regionNames)
    print(paste('rotate',i))
}

for (i in 1:100){

    genes = heatGeneOut(paste0(rotationOut,'/',i,'/Relax/'), heatGenes,2,T)
    heatUp(paste0(outFolder,'/',finalExp),
           paste0(rotationOut,'/',i,'/rotateElim'),
           paste0(rotationOut,'/',i,'/heatmapElim.png'),
           heatProps,
           heatColors,
           heatPalette,
           genes,
           2000,
           2000
           )

    heatUp(paste0(outFolder,'/',finalExp),
           paste0(rotationOut,'/',i,'/rotateForward'),
           paste0(rotationOut,'/',i,'/heatmapForward.png'),
           heatProps,
           heatColors,
           heatPalette,
           genes,
           2000,
           2000
           )
    print(paste('heatmap',i))
    
}



source('rotateCheck.R')
rotateCheck(rotationOut)

source('rotationSelect.R')
rotationSelect(rotationOut, geneOut, rotSelOut)





source('humanBipol.R')
humanBipol(paste0(geneOut,'/Relax/Cortex_CellType'), bipolLoc, bipolOut, paste0(outFolder,'/',finalExp))

# source('coexpAna.R')
# coexpAna(humanDes, genesOut,regionMapping, humanDat,  coexpHOut)


# calculates specificity index as describe in
# https://www.landesbioscience.com/journals/systemsbiomedicine/article/25630/
# takes too much time!
# source('specIndex.R')
# specIndex(paste0(outFolder,'/meltedDesign'),
#          paste0(outFolder,'/finalExp'),
#          geneOutIndex,
#          groupNames)
