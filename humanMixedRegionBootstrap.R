require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
source('readHumanCel.R')

#humanMixedRegionBootstrap = function(humanSoftLoc,hCelDir, hBootOut)

#humanSoft  = read.table('Data/hugeHumanSoft.tsv', sep = '\t', header = T,quote = '')
hBootOut = 'Data/humanBootstrap'
hCelDir = 'humanRegionCel'

humanRegions = read.table('Data/hugeHumanSoft.tsv', header = T, sep = '\t', stringsAsFactors = F)

humanRegions = humanRegions[humanRegions$platform=='GPL5175'
                            & !humanRegions$CDeath=='Cancer',]


require(foreach)
require(doMC)
require(parallel)
cores = 4
# so that I wont fry my laptop
if (detectCores()<cores){ cores = detectCores()}
registerDoMC(cores)


library(oligo)
library(pd.huex.1.0.st.v2)

cels = list.celfiles(hCelDir)

dir.create(hBootOut,showWarnings = F, recursive = T)

size = 15
# foreach(i = 2:100) %dopar% {
for (i in 474:500){
    allRegions = unique(humanRegions$Region)
    subsetInd = vector(length = len(allRegions)*size)
    for (j in 1:len(allRegions)){
        subsetInd[(j*size+1-size):(j*size)] = sample(which(humanRegions$Region %in% allRegions[j]),size)
    }
    subsetRegions = humanRegions[subsetInd,]
    readHumanCel(subsetRegions$GSM,paste0('Data/HumanBootstrap/',i),humanDir='humanRegionCel')
    
}





