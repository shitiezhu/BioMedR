
# for batched cor.test

corr_test <- function(data = NULL, gene = NULL, method = c("pearson", "kendall", "spearman")){
  y <- as.numeric(data[gene,])
  corlist <- apply(data, 1, function(x){
    cor.test(as.numeric(x), y, method = method)
  })
  purrr::map_dfr(corlist, broom::tidy, .id = 'mRNA')
}
