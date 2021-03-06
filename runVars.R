require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
# sourceGithub(oganm,toSource,gemmaAnnotate)

# getGemmaAnnotGoogle('GPL339','Data/GPL339Annotation')

write.design = function(x, file){
    write.table(x,file= file, sep = '\t', quote=F, row.names = F)
}

read.design  = function(x){
    read.table(x,header=T,sep='\t',stringsAsFactors=F,quote="")
}

read.exp = function(x){
    read.csv(x,header = T,stringsAsFactors=F)
}

sepExpr = function(allDataPre){
    for (i in 1:ncol(allDataPre)){
        if ('double'==typeof(allDataPre[,i])){
            expBound = i
            break
        }
    }
    geneData = allDataPre[,1:(expBound-1)]
    exprData = allDataPre[,expBound:ncol(allDataPre)]
    return(list(geneData,exprData))
}

# for loading and normalization -----
desFile='Data/Design.tsv'
xls = F
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'Normalize'
namingCol2='Normalize2.0'
tinyChip = 'mouse430a2.db'
skipNorm = T
# for gene selection -----
geneOut = 'Data/Fold'
geneOutIndex = 'Data/Index'
groupNames = c('GabaDeep','PyramidalDeep','JustPyra')
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
heatProps = c('GabaDeep',
              'MajorType',
              'Reference',
              'Age')
heatGenes = 'GabaDeep'
heatPalette = colorRampPalette(c("darkred",'white', "blue"))(n = 1000)
# long coloring var for heatmap. modify to your heart's content ----
# extra cell types will not cause problems if in list
namingColors = c(method = 'direct',
                 legendLoc = 'bottomleft',
                 cex = 2.3,
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
                 # Gaba = 'firebrick4',
                 Astrocyte = 'yellow',
                 GabaPV = 'firebrick2',
                 Stem = 'blue' ,
                 Ependymal = 'orange',
                 Serotonergic = 'darkolivegreen',
                 Hypocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0',
                 Th_positive_LC = 'blueviolet',
                 GabaVIPReln = 'firebrick4',
                 GabaRelnCalb = 'firebrick3',
                 GabaSSTReln = 'firebrick1',
                 GabaReln = 'firebrick')


ageColors = c(method = 'gradFactor', lo = 'blue', hi = 'red',
              legendLoc = 'bottomright',factorOrd = NA,
              cex = 2.3)

typeColors = c(method = 'def', legendLoc = 'none')

refColors = c(method = 'def', legendLoc = 'topleft', cex = 2.3)
# names are not important
heatColors = list(CellType = namingColors, typeColors = typeColors, Reference = refColors, Age = ageColors)

# heatmap for genes
heatColors2 = heatColors[c(1,3)]
heatColors2[[1]]['cex'] = '3'
heatColors2[[2]]['legendLoc'] = "none"
#heatColors2[[4]]['legendLoc'] = "none"

# for bipol ----
bipolLoc = 'Data/GSE5388_GSE12649_for_Ogan.RData'
bipolOut = 'Data/BipolPC'
bipolOutHeat = 'Data/BipolPC/Heatmaps'