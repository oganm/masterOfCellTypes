library(oligo)
library(pd.huex.1.0.st.v2)
source('readHumanCel.R')
softFile = read.design('Data/hugeHumanSoft.tsv')
softFile = softFile[!is.na(softFile$pH) 
                    & softFile$platform=='GPL5175'
                    & !softFile$CDeath=='Cancer',]


regions = unique(softFile$Region)
dir.create('Data/HumanRegionExpr',showWarnings=F)
for (i in regions[2:length(regions)]){
    regionSet = softFile[softFile$Region == i,]
    readHumanCel(regionSet$GSM,paste0('Data/HumanRegionExpr/',i),humanDir='humanRegionCel')
}


# for white matter and cortex, a special to check their differences.
regionSet = softFile[softFile$Region %in% c('frontal cortex','white matter'),]
readHumanCel(regionSet$GSM,paste0('Data/HumanRegionExpr/','cortex-white'),humanDir='humanRegionCel')

regionSet = softFile[softFile$Region %in% c('frontal cortex','cerebellar cortex'),]
readHumanCel(regionSet$GSM,paste0('Data/HumanRegionExpr/','cortex-cerebellum'),humanDir='humanRegionCel')