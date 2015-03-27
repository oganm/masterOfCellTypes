#treshold detection for absance presence call in RNA seq. Seperated from
# coexistance test due to slower speed

require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

require(runVars.R)

# load rna seq data--------------

rnaSeq = read.table('Data/RNASeq/expression_mRNA_17-Aug-2014.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))

require(foreach)
require(doMC)
require(parallel)
require(mixtools)

cores =16
if (detectCores()<cores){ 
    cores = detectCores()
    print('max cores exceeded')
    print(paste('set core no to',cores))
}
registerDoMC(cores)
# to every gene, fit two gaussians, lower one is counted as non expressed,
# higher is expressed
gaus = foreach(x = 1:nrow(rnaExpAll)) %dopar% {
    print(x)
    tryCatch({
        normalmixEM(rnaExp[x,],maxrestarts=5, verb = F)
    }, error = function(e){
        return(NA)
    })
}
# old non parallel version
# gaus = apply(rnaExp[,],1,function(x){
#     tryCatch({
#         normalmixEM(x,maxrestarts=10)
#     }, error = function(e){
#         return(NA)
#     })
# })
 
# find the probability of each count. set the treshold to the point where the 
# probability of a being in the expressed group is higher than probability of 
# being in the non expressed group.
tresholds = sapply(gaus,function(x){
     if (is.na(x)){
         return(1)
     }
    prob1=pnorm(0:max(x$x),mean = x$mu[1], x$sigma[1],lower.tail=F)
    prob2=pnorm(0:max(x$x),mean = x$mu[2], x$sigma[2],lower.tail=T)
    return(min(which(prob2>prob1)))
 })
 
write.table(tresholds,file='Data/RNASeq/tresholds')

