
# simulated coexpression plot ---------
setwd('Documents/Committee Presentation 1')
require(ggplot2)
red = c(3,3,2)
blue= c(4,2,5)
yellow= c(2,4,2)


fakeGenes = data.frame(1:3,jitter(red),jitter(red),jitter(blue),jitter(blue),jitter(yellow),jitter(yellow))
names(fakeGenes) = c('region','red1','red2','blue1','blue2','yellow1','yellow2')
require(reshape)
meltedGenes = melt(fakeGenes,measure.vars = c('red1','red2','blue1','blue2','yellow1','yellow2'))

manualColor =  scale_colour_manual(name='prop', values = c(red1=toupper('#df5a49'),red2=toupper('#df5a49'),
                                                           blue1 = toupper('#334d5c'), blue2= toupper('#334d5c'),
                                                           yellow1= toupper('#efc94c'), yellow2=toupper('#efc94c')))
svg('4.coexpression.svg')
(ggplot(meltedGenes, aes(x=region,y=value,color=variable,group=variable))+
     geom_point(size = 3)+
     geom_line(size = 2)+
     manualColor+
     scale_x_discrete(breaks = 1:3, labels=c('Region1','Region2','Region3'),name='')+
     theme(legend.position="none")+
     scale_y_continuous(limits=c(1,5.5),name='Relative expressions of genes' , breaks= NULL)+
     theme(axis.title.y = element_text(size=18))+
     theme(axis.text.x = element_text(size = 18))
)

dev.off()

# simulated expression change plot before -------------------

fakeGenes = data.frame(gene = factor(paste0('Gene',1:7)),
                       expression =  runif(7,6,16) , 
                       color = c('blue' , rep('red',3),rep('blue', 3)))

manualColor =   scale_fill_manual(name='prop',
                                  values = c(red = '#DF5A49', blue ='#334D5C'))

svg('1.exprBefore.svg')
(ggplot(fakeGenes, aes (x=gene, y = expression , fill = color)) + 
     manualColor + 
    geom_bar(stat='identity') + 
     scale_x_discrete(name = ''  )+
    scale_y_continuous(limits = c(0,16), breaks = NULL)+ 
     theme(legend.position="none") + 
     theme(axis.title.y = element_text(size=20)) + 
     theme(axis.text.x = element_text(size = 18))     
 )

dev.off()

# simulated expression change after -------------------
svg('2.exprAfter.svg')

fakeGenes$expression[fakeGenes$color=='red'] = fakeGenes$expression[fakeGenes$color=='red'] / 4
(ggplot(fakeGenes, aes (x=gene, y = expression , fill = color)) + 
     manualColor + 
     geom_bar(stat='identity') + 
     scale_x_discrete(name = ''  )+
     scale_y_continuous(limits = c(0,16), breaks = NULL)+ 
     theme(legend.position="none") + 
     theme(axis.title.y = element_text(size=20)) + 
     theme(axis.text.x = element_text(size = 18))     
)

dev.off()

# human data degredation plot ------
# code looks weird because it's taken from qc file. that whole list business
# is not necessary
require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
parent = getParent(2)
source('runVars.R')

setwd(parent)



humanDesign = read.table('../../Data/humanCellTypeDesign.tsv',header=T,sep='\t', stringsAsFactors=F)


gsms = humanDesign[,1]

gsms = strsplit(gsms,',')


platforms = unique(humanDesign$Platform)

require(affy)

affies = lapply( rep("AffyBatch", len(platforms)), new )
celsInFolder = list.files('../../humanCellTypeCel')
celsNoExtension = gsub('[.](C|c)(E|e)(L|l)','',celsInFolder)

for (i in 1){
    toRead = unlist(gsms[humanDesign$Platform %in% platforms[i]])
    relevant = celsInFolder[celsNoExtension %in% toRead]
    affies[i] = ReadAffy(filenames = paste0('../../humanCellTypeCel/',relevant))
}

names(affies) = platforms



setwd('../..')
svg('5.exprAfter.svg')

plotAffyRNAdeg(affies[[1]])
dev.off()

# not working markers plot ----------------
load(bipolLoc)
bpCntScz = aned_high_GSE12649
bpCnt = aned_high_GSE5388
bpCntSczDes = Samples_GSE12649
bpCntDes = Samples_GSE5388
rm(aned_high_GSE12649)
rm(aned_high_GSE5388)
rm(Samples_GSE12649)
rm(Samples_GSE5388)

windowSize = 4
frame = data.frame(t(bpCnt[c('NOV','FAM131A'),4:ncol(bpCnt)]),
           group = gsub('_.*','',names(bpCnt[4:ncol(bpCnt)])))

frame = melt(frame)
svg("6. not working gene expression")
(ggplot(frame, aes(x= variable, y = value))
 + geom_point()  
 + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20),
         axis.title.y = element_text(vjust=0.5, size=20),
         axis.title.x = element_text(vjust=0.5, size=0) ,
         title = element_text(vjust=0.5, size=20),
         strip.text.x = element_text(size = 20),
         axis.text.y = element_text(size = 16))
 +  scale_y_continuous(name="log2 Expression"))

def.off()

