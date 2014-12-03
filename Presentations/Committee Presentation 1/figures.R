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

(ggplot(meltedGenes, aes(x=region,y=value,color=variable,group=variable))+
     geom_point(size = 3)+
     geom_line(size = 2)+
     manualColor+
     scale_x_discrete(breaks = 1:3, labels=c('Region1','Region2','Region3'),name='')+
     theme(legend.position="none")+
     scale_y_continuous(limits=c(1,5.5),name='Relative expressions of genes')+
     theme(axis.title.y = element_text(size=16))
)




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