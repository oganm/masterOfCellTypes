# Reads CEL files from a given directory according to a design file
#essential parts of the file are
#1. the first collumn must include the GSMs
#2. platform must be written to the Platform collumn

require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

sourceGithub(oganm,toSource,gemmaAnnotate)
sourceGithub(oganm,toSource,mergeChips)


readDesignMergeCel = function (desFile, namingCol, celRegex, celDir,tinyChip, outFolder){
    #always have gsms in the first collumn

   
    #eval( expr = parse( text = getURL(
    #    "https://raw.githubusercontent.com/oganm/toSource/master/mergeChips.R",
    #    ssl.verifypeer=FALSE) ))
    #only works in windows
    #source('https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r')
    #source('https://raw.githubusercontent.com/oganm/toSource/master/mergeChips.R')
    require(affy)
    require(compare)



    design = read.table(desFile,quote='',header=T,sep='\t')

    design = design[!is.na(design[,namingCol]),]

    gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))
    if (any(lapply(gsms,len)==0)){
        print('No cel names could be captured from the following rows')
        print((1:len(gsms))[lapply(gsms,len)==0])
    }

    platforms = unique(design$Platform)

    affies = lapply( rep("AffyBatch", len(platforms)), new )
    names(affies) = platforms

    #create a vector of affy objects based on how many chips are there just works with 2 right now
    #implementation of a third generation is straightforward. in fact i could have done it instead
    #of writing this but oh well...
    for (i in 1:len(affies)){
        celsInFolder = list.files(
            paste0(celDir,'/',platforms[i]))
        celsNoExtension = gsub('[.](C|c)(E|e)(L|l)','',celsInFolder)
        gsms = regmatches(design[design$Platform == platforms[i], 1], gregexpr(celRegex, design[design$Platform == platforms[i], 1],perl=T))
        gsms = unlist(gsms)

        relevant = celsInFolder[celsNoExtension %in% gsms]

        affies[i] = ReadAffy(filenames = paste0(celDir,'/',platforms[i],'/',relevant))
    }

    newNormalized = mergeChips(affies[[1]],affies[[2]])
    
    aned = gemmaAnnot(newNormalized, 'Data/GPL339Annotation')
    
    aned = aned[!aned$Gene.Symbol == '',]
    
    

    header = gsub('.cel', '', gsub('.CEL','', colnames(aned)[4:ncol(aned)]))
    colnames(aned) = c(colnames(aned)[1:3], header)
    dir.create(outFolder, recursive = T,showWarnings=F)
    write.csv(aned, paste0(outFolder,"/rmaExp.csv"), row.names=FALSE)
    #boxplot(aned[,4:ncol(aned)])
    gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))

    #gsms = regmatches(df[, 1], gregexpr("(GSM\\d\\d\\d\\d\\d(\\d|))|(PC\\d....)|(Y+.*?((?=(,))|\\d+))|(((?<=:)|(?<=,))A\\d.*?30A)|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))", df[, 1],perl=T))


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

    write.table(newDesign, paste0(outFolder,"/meltedDesign.tsv"), row.names=FALSE,sep = '\t', quote=F)
}


# for changes in original design file that only involves naming. do not run the whole thing again
meltDesign = function(desFile, namingCol, celRegex, exprFile, outFile){
    expr = read.csv(exprFile , header = T)
    header =  colnames(expr)[4:ncol(expr)]

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

