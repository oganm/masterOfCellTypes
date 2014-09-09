input = 'Data2/normalizedDesign'
output = input
expr = 'Data2/mostVariableQuantileNormalized'

allDataPre = read.csv(expr, header = T)
design = read.table(input,header=T,sep='\t')

geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]

design = design[match(colnames(exprData),gsub('[+]','.',gsub('-','.',design$sampleName))),]


Xist = which(geneData$Gene.Symbol %in% 'Xist')

sex = rep('varys', ncol(exprData))

sex[exprData[Xist,]>=7] = 'female'
sex[exprData[Xist,]<7] = 'male'


newDesign = cbind(design,sex)
newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

write.table(newDesign,file=output, quote=F, sep= '\t',row.names=F) 
