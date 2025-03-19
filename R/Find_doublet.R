Find_doublet <- function(data){

  data = JoinLayers(data)
  sweep.res.list <- DoubletFinder::paramSweep(data, PCs = 1:20, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- DoubletFinder::doubletFinder(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data = subset(data, subset = doublet_info == "Singlet")
  return(data)
}
