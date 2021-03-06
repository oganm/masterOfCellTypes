Quality Control Report
========================================================





This file contains the quality control information on the chips we used and choose not to use. 

We were able to find limited amounts of human cell type specific data but they were shown to be of very poor quality. The sources and GEO IDs of these datasets can be found from [this](/Data/humanCellTypeDesign.tsv) file

Human Cell Type Data
-------------
Human data was found on 3 chips `platforms r`. Below, the percentages of MM>PM signals on each chip group can be found

```r
lapply(affies,function(x){
    mean(mm(x)>pm(x))
})
```

```
## $`Affymetrix U133 Plus 2.0`
## [1] 0.3319178
## 
## $`Affymetrix U133_X3P`
## [1] 0.3612196
## 
## $`Affymetrix U133A`
## [1] 0.4254485
```

Per chip ratio of MM>PM vs overall pm+mm intensity


```r
lapply(affies,function(x){
data = data.frame(MMPM = rep(0,len(x)), intensity = rep(0,len(x)))
for (i in 1:len(x)){
    data$MMPM[i] = mean(mm(x[,i]) > pm(x[,i]))
    data$intensity[i] = mean(pm(x[,i])+mm(x[,i]))
}
plot(data$MMPM,data$intensity)

}
)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png) 

```
## $`Affymetrix U133 Plus 2.0`
## NULL
## 
## $`Affymetrix U133_X3P`
## NULL
## 
## $`Affymetrix U133A`
## NULL
```


A quite messy plots of density distributions. Note that there are multiple batches in the same chip 

```r
ground = lapply(affies,hist)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) ![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) ![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png) 

And RNA degredation plots with mean slopes

```r
ground = lapply(affies, function(x){
    deg = AffyRNAdeg(x)
    print(mean(deg$slope))
    plotAffyRNAdeg(deg)
})
```

```
## [1] 8.691982
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```
## [1] 6.182243
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

```
## [1] 3.976004
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png) 

In the first plot we notice a single sample with quite a low slope and decided to see the raw image file for that chip. Which turned out to be an unusually homogenous image


```r
deg = AffyRNAdeg(affies[[1]])
image(affies[[1]][,order(deg$slope)[1]])
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 


A similar outlier also exists in the 3rd chip type

```r
deg = AffyRNAdeg(affies[[3]])
image(affies[[3]][,order(deg$slope)[1]])
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 


Mouse Cell Type
------------------
Lets proceed to our lovely mouse data, that has two platforms






```r
lapply(affies,function(x){
    mean(mm(x)>pm(x))
})
```

```
## $`Affymetrix 430 2.0`
## [1] 0.2769139
## 
## $`Affymetrix MOE430A`
## [1] 0.2916909
```

```r
lapply(affies,function(x){
data = data.frame(MMPM = rep(0,len(x)), intensity = rep(0,len(x)))
for (i in 1:len(x)){
    data$MMPM[i] = mean(mm(x[,i]) > pm(x[,i]))
    data$intensity[i] = mean(pm(x[,i])+mm(x[,i]))
}
plot(data$MMPM,data$intensity)

}
)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-2.png) 

```
## $`Affymetrix 430 2.0`
## NULL
## 
## $`Affymetrix MOE430A`
## NULL
```

```r
ground = lapply(affies,hist)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-3.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-4.png) 

```r
ground = lapply(affies, function(x){
    deg = AffyRNAdeg(x)
    print(mean(deg$slope))
    plotAffyRNAdeg(deg)
})
```

```
## [1] 5.437618
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-5.png) 

```
## [1] 4.84261
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-6.png) 

