require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

read.design  = function(x){
    read.table(x,header=T,sep='\t')
}

read.exp = function(x){
    read.csv(x,header = T)}

# for loading and normalization -----
desFile='Data/Design.tsv'
xls = F
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data'
namingCol = 'Normalize'
tinyChip = 'mouse430a2.db'
skipNorm = T
# for gene selection -----
geneOut = 'Data/Fold'
geneOutIndex = 'Data/Index'
groupNames = c('CellType', 'GabaDeep', 'forContanim')
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
                 Hypocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0',
                 Th_positive_LC = 'blueviolet')


ageColors = c(method = 'gradFactor', lo = 'blue', hi = 'red',
              legendLoc = 'bottomright',factorOrd = NA,
              cex = 1.5)

typeColors = c(method = 'def', legendLoc = 'none')

refColors = c(method = 'def', legendLoc = 'topleft', cex = 1.5)
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