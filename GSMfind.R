GSMfind = function(GSE, regex=''){
    require(RCurl)
    page = getURL(paste0('www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', GSE))
    gsms = regmatches(page,gregexpr(paste0('GSM[0-9]*?(?=<.*\n.*?',regex,'</td)'),page,perl=T))
    return(gsms)
}
