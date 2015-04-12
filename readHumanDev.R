library(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
sourceGithub(OganM,toSource,gemmaAnnotate)

library(oligo)
library(stringr)
meta = read.table('BrainDev/Metadata_GSE25219', sep='\t', stringsAsFactors=F)

regionNames = read.table('BrainDev/Region_codes', stringsAsFactors=F)

series = list(OFC = 'OFC',
              DFC = 'DFC',
              VFC = 'VFC', 
              MFC = 'MFC',
              M1C = 'M1C',
              OFC_full = c('FC','OFC'),
              DFC_full = c('DFC','DFC'),
              VFC_full = c('FC','VFC'),
              MFC_full = c('FC','MFC'),
              M1C_full = c('FC','M1C'),
              S1C = 'S1C',
              IPC = 'IPC',
              S1C_full = c('IPC','S1C'),
              IPC_full = c('IPC','IPC'),
              A1C = 'A1C',
              STC = 'STC',
              ITC = 'ITC',
              A1C_full = c('TC','A1C'),
              STC_full = c('TC','STC'),
              ITC_full = c('TC','ITC'),
              V1C = 'V1C',
              V1C_full = c('OC','V1C'),
              HIP = 'HIP',
              AMY = 'AMY',
              STR = 'STR',
              MD = 'MD',
              MD_full = c('DIE', 'MD'),
              CBC = 'CBC',
              CBD_full = c('URL', 'CBC')
              )

cels = celFiles('BrainDev/',listGzipped=T)
              
dir.create('Data/DevelopData')
for (i in 3:len(series)){
    relevant  = rn(meta)[meta$region %in% series[[i]]]
    toRead = cels[str_extract(cels,perl('GSM.*?(?=_)')) %in% relevant]
    oligoRaw = oligo::read.celfiles(paste0('BrainDev/',toRead))
    exonTS = oligo::rma(oligoRaw , target='core')
    featureData(exonTS) <- getNetAffx(exonTS, "transcript")
    aned = gemmaAnnotOligo(normalized=exonTS,chipFile='Data/GemmaAnnots/GPL5175')
    colnames(aned)[grepl('_',cn(aned))]= str_extract(cn(aned)[grepl('_',cn(aned))], perl('^.*?(?=_)'))
    write.csv(aned,paste0('Data/DevelopData/',names(series)[i]),row.names=F)
}
              
              