require('biomaRt')
humanMouseOrthologue = function(humanGenes){
    # human gene symbol in, human symbol with orthologue mouse symbol out
    humanMart = useMart("ensembl",dataset="hsapiens_gene_ensembl")

    mouseMart = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

    symbolID = getBM(c('hgnc_symbol','ensembl_gene_id'),
                     filters = 'hgnc_symbol', values = humanGene$Gene_Symbol,
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