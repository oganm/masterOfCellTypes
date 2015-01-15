require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

sourceGithub(user=OganM,repo= toSource, script=homologene)

require(reshape2)
require(ggplot2)
# it takes in a list of estimates or a single vector of estimates. and plots things
plotEstimates = function(estimates,groups,plotNames, sigTest =  wilcox.test,
                         pAdjMethod = p.adjust.methods,
                         correction = c('columnar', 'all'),
                         comparisons = 'all',
                         sigTreshold = 0.05){
    groupNames = unique(groups)
    
    if (typeof(estimates)!='list'){
        estimates = list(estimates)
    }
    
    if (comparisons == 'all'){
        comparisons = combn(groupNames,2)
    }
    
    # create p value lists for correction
    pList = matrix(data = NA, ncol = ncol(comparisons), nrow = len(estimates))
    for (i in 1:len(estimates)) {
        for (j in 1:ncol(comparisons)){
            pList[i, j] = sig.test(estimates[[i]][groups %in% comparisons[1,j]],
                                   estimates[[i]][groups %in% comparisons[2,j]])$p.value
        }
    }
    # p value adjustment
    if (correction[1] == 'columnar'){
        for (i in 1:ncol(pList)){
            pList[,i] = p.adjust(pList[,i], pAdjMethod)
        }
    }
    
    if (correction[1]=='all'){
        pList = matrix(p.adjust(cbind(pList,pList)),nrow = len(estimates))
    }
    
    
    # plotting of things
    for (i in 1:len(estimates)){
        frame = data.frame(PC1 = estimates[[i]], group = groups)
        windowUp = max((frame$PC1)) + 1
        windowDown = min((frame$PC1)) - 0.5
        
        # prepare significance text
        sigText = apply(comparisons, 2, paste0, collapse = '/')
        
        sigText = paste0(sigText,': ',sprintf('%.5f',pList[i,]),collapse = '\n')
        
        lePlot = ggplot(frame,aes(x=group, y = PC1)) +
            geom_violin( color="#C4C4C4", fill="#C4C4C4") + 
            geom_boxplot(width=0.1,fill = 'lightblue') + 
            geom_point(size = 3) + 
            ggtitle(names(estimates)[i]) + 
            scale_y_continuous(limits=c(windowDown, windowUp),
                               name="Relative estimate of cell type amounts") +
            theme_bw() + 
            theme(axis.text.x  = element_text(size=25),
                  axis.title.y = element_text(vjust=0.5, size=25),
                  axis.title.x = element_text(vjust=0.5, size=0) ,
                  title = element_text(vjust=0.5, size=25),
                  axis.text.y = element_text(size = 13)) + 
            annotate('text', x = 0.1 , y = windowUp ,
                     label = sigText ,
                     hjust = 0, vjust=1, size = 4.5)
        (lePlot)
        ggsave(plotNames[i])
        
    }
    frame = data.frame(PC1 = estimates, group = groups)
    groupNames = unique(groups)
    
}



# frame = data.frame(PC1 = pca$x[,1], group = groups)
# windowSize = max(abs(frame$PC1)) + 0.5
# groupNames = unique(groups)


# estimates relative abundances of cell types based on the PC1s of gene expressions accross samples
# controlBased is the name of the control group in 'groups'. if given, PC calculations will be
# done on controls and rotations found will be applied to other samples.
# groups variable should be provided if the plot is groupBased and/or if 
# controlBased is set to a name of a group
cellTypeEstimate = function(exprData, 
                            genes,
                            geneColName = 'Gene.Symbol',
                            geneTransform = function(x){mouse2human(x)$humanGene}, 
                            groups= NA,
                            controlBased = NA, 
                            tableOut=NA,
                            indivGenePlot = NA,
                            plotType = c('groupBased','cummulative')){
    
    if (typeof(genes)!='list'){
        genes = list(genes)
    }
    
    output = vector(mode = 'list', length=len(genes))
    for (i in 1:len(genes)){
        if(!is.na(geneTransform)){
            genes[[i]] = geneTransform(genes[[i]])
        }
        
        relevantData = exprData[exprData[, geneColName] %in% genes[[i]],]
        rownames(relevantData) = relevantData[, geneColName]
        relevantExpr = relevantData[4:len(relevantData)]
        
        
        if (!is.na(indivGenePlot[1])){
            indivGenes = data.frame(melt(relevantExpr)[,2],rownames(relevantExpr))
            names(indivGenes) = c('expression', 'gene')
            
            switch(plotType[1],
                   groupBased = {
                       indivGenes =  data.frame(indivGenes,group = groups)
                   },
                   cummulative = {
                       indivGenes =  data.frame(indivGenes,group = '')
                   })
            
            p = ggplot(indivGenes,aes(y = expression, x = group )) +
                facet_wrap('gene') + 
                geom_point() +
                theme(axis.text.x  = element_text( size=20),
                      axis.title.y = element_text(vjust=0.5, size=20),
                      axis.title.x = element_text(vjust=0.5, size=0) , 
                      title = element_text(vjust=0.5, size=20))+
                scale_x_discrete(name = '')+
                scale_y_discrete(name = 'log2 Expression')
            (p)
            ggsave(filename = indivGenePlot[i])
        }
        
        
        if (is.na(controlBased)){
            pca = prcomp(t(relevantExpr), scale = T)
        } else{
            pca = prcomp(t(relevantExpr[groups %in% controlBased]), scale = T)
        }
        
        pca$rotation = pca$rotation * ((sum(pca$rotation[,1])<0)*(-2)+1)
        
        
        
        if (!is.na(tableOut[1])){
            
            # add variation explained as a comment on top
            file.create(tableOut[i])
            fileConn = file(tableOut[i])
            writeLines(paste0('# Variation explained: ', 
                              paste(summary(pca)$importance[2,], collapse=' ')), 
                       fileConn)
            close(fileConn)
            write.table(pca$rotation[,1,drop=F],
                        file = tableOut[i],
                        quote = F, row.names = T, col.names = F, sep='\t',
                        append = T)
        }
        
        
        pca$x = t(as.matrix(t(scale(t(relevantExpr))))) %*% as.matrix(pca$rotation)
        
        output[[i]]=(pca$x[,1])
    }
    names(output) = names(genes)
    return(output)
}