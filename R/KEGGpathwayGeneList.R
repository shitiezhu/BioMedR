KEGGpathwayGeneList <- function(dbentries = NULL){
  require(magrittr)
  genenames <- KEGGREST::keggGet(dbentries = dbentries)[[1]]$GENE %>%
    split(.,rep(1:2, length(.)/2)) %>%
    cbind.data.frame() %>%
    setNames(c("ENTREZID","Gene")) %>%
    dplyr::separate(Gene, c("SYMBOL", "Names"), sep = ";")
  genenames
}
