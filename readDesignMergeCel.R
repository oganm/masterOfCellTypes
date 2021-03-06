# Reads CEL files from a given directory according to a design file
#essential parts of the file are
#1. the first collumn must include the GSMs
#2. platform must be written to the Platform collumn
require(stringr)
require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

sourceGithub(oganm,toSource,gemmaAnnotate)
sourceGithub(oganm,toSource,mergeChips)


readMouseCel = function(GSMs,mouseDir='cel',file=NA){
    require(affy)
    require(compare)
    allCells = affy::list.celfiles(mouseDir, recursive = T)
    files = sapply(GSMs,function(x){grep(paste0(x,'[.](C|c)(E|e)(L|l)$'),allCells,value=T)})
    platforms = unique(dirname(unlist(files)))
    affies =  lapply( rep("AffyBatch", len(platforms)), new )
    
    for (i in 1:len(affies)){
        affies[i] = ReadAffy(filenames = paste0('cel/',grep(platforms[i],files, value=T)))
    }
    
    if (len(affies)>1){newNormalized = mergeChips(affies[[1]],affies[[2]])
    } else {newNormalized = rma(affies[[1]])}
    
    if('GPL339' %in% platforms){
        aned = gemmaAnnot(newNormalized, 'Data/GemmaAnnots/GPL339')
    } else{
        aned = gemmaAnnot(newNormalized, 'Data/GemmaAnnots/GPL1261')
    }
    aned = aned[!aned$Gene.Symbol == '',]
    
    names(aned) = gsub('[.](C|c)(E|e)(L|l)','',names(aned))
    
    if (!is.na(file)){
        write.csv(aned, file, row.names=FALSE)
    }
    invisible(aned)
    
}

readDesignMergeCel = function (desFile, namingCol, celRegex, celDir,tinyChip, expFile,desOut){
    #always have gsms in the first collumn

    require(affy)
    require(compare)

    design = read.table(desFile,quote='',header=T,sep='\t')

    design = design[!is.na(design[,namingCol]),]

    gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))
    if (any(lapply(gsms,len)==0)){
        print('No cel names could be captured from the following rows')
        print((1:len(gsms))[lapply(gsms,len)==0])
    }

    aned = readMouseCel(unlist(gsms),celDir)
        
    dir.create(outFolder, recursive = T,showWarnings=F)
    
    list[genes,exp] = sepExpr(aned)
    #boxplot(aned[,4:ncol(aned)])
    gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))

    #gsms = regmatches(df[, 1], gregexpr("(GSM\\d\\d\\d\\d\\d(\\d|))|(PC\\d....)|(Y+.*?((?=(,))|\\d+))|(((?<=:)|(?<=,))A\\d.*?30A)|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))", df[, 1],perl=T))

    header = names(exp)

    indexes = vector()
    for (i in 1:length(header)){
        indexes = c(indexes, findInList(header[i], gsms))
    }
    header[!header %in% unlist(gsms)]
    newDesign = data.frame(sampleName = header, originalIndex = indexes, design[indexes,])
    colnames(newDesign) = c('sampleName','originalIndex',names(design))
    #newDesign$originalIndex = as.numeric(newDesign$originalIndex)

    #newDesign$age = gsub('~', '', newDesign$age)
    #newDesign$age = gsub('P', '', newDesign$age)
    #newDesign$age = gsub('7-8', '7.5', newDesign$age)
    #newDesign$age[grepl('(precise age not given)',newDesign$age)] = 60
    #newDesign$age = as.numeric(newDesign$age)
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

    write.table(newDesign, desOut, row.names=FALSE,sep = '\t', quote=F)
    
    #list[genes,expr]=sepExpr(aned)
    
    #expr = expr[match(make.names(design$sampleName)]
    
    exp = exp[,match(make.names(newDesign$sampleName),make.names(colnames(exp)))]
    aned = cbind(genes,exp)
    write.csv(aned, expFile, row.names=FALSE)
    
}


# for changes in original design file that only involves naming. do not run the whole thing again
meltDesign = function(desFile, namingCol, celRegex, exprFile, outFile){
    expr = read.csv(exprFile , header = T)
    list[gene,expres] = sepExpr(expr)
    header =  colnames(expres)

    design = read.table(desFile,quote='',header=T,sep='\t')
    design = design[!is.na(design[,namingCol]),]
    gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))
    if (any(lapply(gsms,len)==0)){
        print('No cel names could be captured from the following rows')
        print((1:len(gsms))[lapply(gsms,len)==0])
    }

    indexes = vector()
    for (i in 1:length(header)){
        indexes = c(indexes, findInList(header[i], lapply(gsms,make.names)))
    }


    if (len(header[!header %in% make.names(unlist(gsms))])>0){
        print("what the fuck man! I can't find some GSMs in your design file")
        print(header[!header %in% make.names(unlist(gsms))])
    }

    newDesign = data.frame(sampleName = header, originalIndex = indexes, design[indexes,])
    colnames(newDesign) = c('sampleName','originalIndex',names(design))
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]
    write.table(newDesign, paste0(outFile), row.names=FALSE,sep = '\t', quote=F)

}

