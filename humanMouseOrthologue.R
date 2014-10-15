require('biomaRt')
humanMouseOrthologue = function(humanGenes){
    # human gene symbol in, human symbol with orthologue mouse symbol out
    humanMart = useMart("ensembl",dataset="hsapiens_gene_ensembl")

    mouseMart = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

    symbolID = getBM(c('hgnc_symbol','ensembl_gene_id'),
                     filters = 'hgnc_symbol', values = humanGenes,
                     mart = humanMart)

    ortho = getBM(c('mmusculus_homolog_ensembl_gene','ensembl_gene_id'),
                  filters = 'ensembl_gene_id', values = symbolID$ensembl_gene_id,
                  mart = humanMart)

    mouseSymbol = getBM(c('ensembl_gene_id','mgi_symbol'),
                        filters = 'ensembl_gene_id', values = ortho$mmusculus_homolog_ensembl_gene,
                        mart = mouseMart)
    names(mouseSymbol) = c('mmusculus_homolog_ensembl_gene','mgi_symbol')
    fullOrtho = merge(merge(symbolID,ortho),mouseSymbol)
    fullOrtho = fullOrtho[,3:4]
    return(fullOrtho)
}

mouseHumanOrthologue = function(mouseGenes){
    humanMart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    mouseMart = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
    
    symbolID = getBM(c('mgi_symbol', 'ensembl_gene_id'),
                     filters = 'mgi_symbol', values = mouseGenes,
                     mart = mouseMart)
    
    ortho = getBM(c('hsapiens_homolog_ensembl_gene','ensembl_gene_id'),
                  filters = 'ensembl_gene_id', values = symbolID$ensembl_gene_id,
                  mart = mouseMart)
    humanSymbol = getBM(c('ensembl_gene_id', 'hgnc_symbol'),
                        filters = 'ensembl_gene_id',values = ortho$hsapiens_homolog_ensembl_gene,
                        mart = humanMart)
    names(humanSymbol) = c('hsapiens_homolog_ensembl_gene' , 'hgnc_symbol')
    fullOrtho = merge(merge(symbolID,ortho),humanSymbol)
    fullOrtho = fullOrtho[,3:4]
    return(fullOrtho)
}