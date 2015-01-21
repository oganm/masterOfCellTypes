pH correlation report
========================================================


This file includes an analysis of correlations of gene expression in brain to pH in different regions of the brain


Loading of the data and preallocation in a probe-based manner


Gettind the correlations


Looking at the correlations of upper and lower %99 interval

```r
print('higher')
```

```
## [1] "higher"
```

```r
quant = quantile(abs(corrs),0.99)
keep = apply(abs(corrs),1,function(x){any(x>quant)})
p=heatmap.2(corrs[keep,], trace = "none",#breaks = seq(-1,1,0.01),
            col = heatPaletteWideWhite, cexCol=1,dendrogram = 'column',symbreaks=T)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```r
quant = quantile(abs(corrs),0.01)
keep = apply(abs(corrs),1,function(x){any(x>quant)})
print('lower')
```

```
## [1] "lower"
```

```r
p=heatmap.2(corrs[keep,], trace = "none",#breaks = seq(-1,1,0.01),
            col = heatPaletteWideWhite, cexCol=1,dendrogram = 'column',symbreaks=T)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) 


enrichment analysis for lilah's ![extreme](http://thumbs.dreamstime.com/x/extreme-adventure-mud-logo-graphic-text-20083928.jpg) genes in the upper quantile


```r
quant = quantile(abs(corrs),0.99)
keep = apply(abs(corrs),1,function(x){any(x>quant)})
x=c(sum(genes_GSE538 %in% rownames(corrs[keep,])), 
    sum(genes_GSE12649 %in% rownames(corrs[keep,])))

dhyper(sum(genes_GSE538 %in% rownames(corrs[keep,])),
       sum(rownames(corrs) %in% genes_GSE538),
       nrow(corrs) - sum(rownames(corrs) %in% genes_GSE538),
       sum(keep))
```

```
## [1] 2.525904e-06
```

```r
dhyper(sum(genes_GSE12649 %in% rownames(corrs[keep,])),
       sum(rownames(corrs) %in% genes_GSE12649),
       nrow(corrs) - sum(rownames(corrs) %in% genes_GSE12649),
       sum(keep))
```

```
## [1] 0.1247169
```

Heatmaps of extreme genes

```r
heatmap.2(corrs[rownames(corrs) %in% genes_GSE538,] , trace = "none",
          col = heatPaletteWideWhite, cexCol=1,dendrogram = 'column')
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```r
heatmap.2(corrs[rownames(corrs) %in% genes_GSE12649,] , trace = "none",
          col = heatPaletteWideWhite, cexCol=1,dendrogram = 'column')
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 


pH correlations of all genes accross regions

```r
frame = melt(abs(corrs))
frame$Var2 = factor (frame$Var2,levels = names(sort(apply(abs(corrs),2,mean),decreasing=T)))

lePlot = ggplot(frame,aes(x=Var2, y = value)) +
    geom_violin( color="#C4C4C4", fill="#C4C4C4") + 
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    scale_y_continuous(name="Correlations to pH") +
    theme_bw() + 
    theme(axis.text.x  = element_text(size=20,angle=90),
          axis.title.y = element_text(vjust=0.5, size=25),
          axis.title.x = element_text(vjust=0.5, size=0) ,
          title = element_text(vjust=0.5, size=25),
          axis.text.y = element_text(size = 13))

(lePlot)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 


pH correlations of top %99 quantile


```r
quant = quantile(abs(corrs),0.99)
keep = apply(abs(corrs),1,function(x){any(x>quant)})

frame = melt(abs(corrs[keep,]))
frame$Var2 = factor (frame$Var2,levels = names(sort(apply(abs(corrs),2,mean),decreasing=T)))

lePlot = ggplot(frame,aes(x=Var2, y = value)) +
    geom_violin( color="#C4C4C4", fill="#C4C4C4") + 
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    scale_y_continuous(name="Correlations to pH") +
    theme_bw() + 
    theme(axis.text.x  = element_text(size=20,angle=90),
          axis.title.y = element_text(vjust=0.5, size=25),
          axis.title.x = element_text(vjust=0.5, size=0) ,
          title = element_text(vjust=0.5, size=25),
          axis.text.y = element_text(size = 13))

(lePlot)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 
