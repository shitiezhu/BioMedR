
# countdata <- as.matrix(openxlsx::read.xlsx("gene_count_matrix - HCT116 Oxa vs Mock.xlsx", rowNames = TRUE))
# group <- as.factor(c("OX", "OX", "NC", "NC"))
#
# design <- model.matrix(~0+group)
# colnames(design) <- gsub("group", "", colnames(design))
# design
#
# contr.matrix <- makeContrasts(
#   OXvsNC = OX - NC,
#   levels = colnames(design))
#
# if(!dir.exists("results")) {dir.create("results")}
#
# TZ_DEG(countdata = countdata, group = group,
#        contr.matrix = contr.matrix, design = design,
#        cut.fc = 2, cut.pval = 0.05, bartopn = 10, enrich = FALSE)

TZ_DEG <- function(countdata = countdata, group = group,
                   contr.matrix = contr.matrix, design = design,
                   cut.fc = 2, cut.pval = 0.05,bartopn = 10, enrich = FALSE) {

  library(limma)
  library(Glimma)
  library(edgeR)
  # BiocManager::install("Homo.sapiens")
  library(Homo.sapiens)

  DGEdata <- edgeR::DGEList(countdata)
  dim(DGEdata)
  colnames(DGEdata)

  DGEdata$samples$group <- group


  lcpm <- cpm(DGEdata, log=TRUE)

  table(rowSums(DGEdata$counts==0)==4)


  keep.exprs <- filterByExpr(DGEdata, group=group)
  DGEdata_keep <- DGEdata[keep.exprs, keep.lib.sizes=FALSE]
  dim(DGEdata_keep)

  DGEdata_keep <- calcNormFactors(DGEdata_keep, method = "TMM")
  DGEdata_keep$samples$norm.factors



  lcpm_keep <- cpm(DGEdata_keep, log=TRUE)

  col.group <- group
  levels(col.group) <-RColorBrewer::brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)

  MDSplot <- plotMDS(lcpm_keep, labels=group, col = col.group)
  MDSplot

  design = design

  contr.matrix = contr.matrix

  v <- voom(DGEdata_keep, design, plot=FALSE)
  v

  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  # plotSA(efit, main="Final model: Mean-variance trend")

  summary(decideTests(efit))

  library(tidyverse)

  geneid <- rownames(DGEdata_keep)
  genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL"),
                                 keytype="ENSEMBL") %>% na.omit() %>%
    filter(!duplicated(SYMBOL))

  counts_for_GSEA <- DGEdata_keep$counts %>% as.data.frame() %>%
    rownames_to_column("ENSEMBL") %>%
    inner_join(genes, by = "ENSEMBL") %>%
    select(-1) %>%
    rename (sample=SYMBOL ) %>%
    mutate(DESCRIPTION = rep("na", nrow(.))) %>%
    .[,c(5,6,1:4)]

  if(!dir.exists("results")) {dir.create("results")}
  write.table(counts_for_GSEA, "results/count matrix.txt",sep = "\t", row.names = F)

  cut.fc = cut.fc
  cut.pval = cut.pval

  DEGresults <- topTreat(efit, coef=1, n=Inf) %>%
    rownames_to_column("ENSEMBL") %>%
    inner_join(genes, by = "ENSEMBL") %>%
    mutate(regulation = case_when(
      logFC > log2(cut.fc) & P.Value < cut.pval ~ "Up_reg",
      logFC < -log2(cut.fc) & P.Value < cut.pval ~ "Down_reg",
      abs(logFC) <= log2(cut.fc) | P.Value >= cut.pval ~ "Not_reg"
    )) %>%
    mutate_at(vars(regulation), as.factor) %>%
    inner_join(as.data.frame(lcpm_keep) %>% rownames_to_column("ENSEMBL"), by = "ENSEMBL")

  write.csv(DEGresults, "results/DEG_analysis_tiezhu.csv", row.names = F)

  table(DEGresults$regulation)

  volcano <- ggplot(DEGresults, aes(logFC, -log10(P.Value), color = regulation))+
    geom_point(size = 2, alpha = 0.6)+
    scale_color_manual(values = c("blue", "grey", "red"))+
    # geom_point(data = filter(OXvsNC, SYMBOL == "PHGDH"), size = 2, color = "black")+
    geom_point(data = top_n(filter(DEGresults, regulation == "Down_reg"), 10, -P.Value), pch = 21, color = "black")+
    geom_point(data = top_n(filter(DEGresults, regulation == "Up_reg"), 10, -P.Value), pch = 21, color = "black")+
    ggrepel::geom_text_repel(data = top_n(filter(DEGresults, regulation == "Up_reg"), 10, -P.Value),
                             aes(label = SYMBOL),max.overlaps = 20)+
    ggrepel::geom_text_repel(data = top_n(filter(DEGresults, regulation == "Down_reg"), 10, -P.Value),
                             aes(label = SYMBOL),max.overlaps = 20)+
    # ggrepel::geom_text_repel(aes(label = ifelse(SYMBOL == "PHGDH", SYMBOL, NA)), color = "black")+
    geom_vline(xintercept = c(-log2(cut.fc), log2(cut.fc)), lty = 2)+
    geom_hline(yintercept = -log10(cut.pval), lty = 2)+
    ggprism::theme_prism(base_fontface = "plain",base_size = 14,base_line_size = 0.5)+
    theme(legend.position = "right")+
    scale_y_continuous(expand = c(0,0), limits = c(0, max(-log10(DEGresults$P.Value) + 1)))+
    labs(x = "log2FC(treatvscontrol)", y = "-log10(p-value)", color = "")
  volcano

  ggsave("results/volcano plot.pdf", volcano, height = 6, width = 7)

  bartopn = bartopn

  p1 <- ggplot(data = top_n(filter(DEGresults, regulation == "Up_reg"), bartopn, -P.Value),
               aes(x = reorder(SYMBOL, logFC), y = logFC))+
    geom_bar(stat='identity',width=0.7, fill = "#A03E29")+
    theme_classic(base_size = 14)+
    theme(axis.text.y = element_text(angle=0, hjust=1,vjust=.5,
                                     color = "black",size = 10))+
    labs(x="", y ="log2FC")+
    # geom_hline(yintercept = log2(2), color = "#00548F", size = .6)+
    coord_flip()

  p2 <- ggplot(data = top_n(filter(DEGresults, regulation == "Down_reg"), bartopn, -P.Value),
               aes(x = reorder(SYMBOL, logFC), y = logFC))+
    geom_bar(stat='identity',width=0.7, fill = "#00548F")+
    theme_classic(base_size = 14)+
    theme(axis.text.y = element_text(angle=0, hjust=1,vjust=.5,
                                     color = "black",size = 10))+
    labs(x="", y ="log2FC")+
    # geom_hline(yintercept = -log2(1.5), color = "#A03E29", size = .6)+
    scale_x_discrete(position = "top")+
    coord_flip()

  top30_up_down <- cowplot::plot_grid(p2, p1, labels = c("a", "b"), align = "h")

  ggsave("results/top_up_down.pdf", top30_up_down, height = 6, width = 7)

  #############################enrichment
  enrich = enrich
  if(enrich == TRUE){
    library (clusterProfiler, quietly = T)
    library(enrichplot, quietly = T)

    genelist_up <- filter(DEGresults, regulation == "Up_reg")$logFC
    names(genelist_up) <- filter(DEGresults, regulation == "Up_reg")$SYMBOL

    GeneList_up <- bitr(names(genelist_up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    head(GeneList_up)

    gene_up <- GeneList_up$ENTREZID

    ego_all_up <- enrichGO(gene = gene_up,
                           OrgDb = "org.Hs.eg.db",
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.5,
                           readable = TRUE)

    ego_bar_up <- barplot(ego_all_up,showCategory=20, split= "ONTOLOGY")+
      facet_grid(ONTOLOGY~.,scale="free")
    # scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    # theme_bw()
    ggsave("results/ego_bar_up.pdf", ego_bar_up, height = 10, width = 10)

    GO_all_up <- as.data.frame(ego_all_up)

    write.csv(GO_all_up, file = "results/GO_all_up.csv")


    kk_up <- enrichKEGG(gene = gene_up,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)
    kk_up <- setReadable(kk_up, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db")

    KEGG_up <- as.data.frame(kk_up)

    write.csv(KEGG_up, file = "results/KEGG_up.csv")

    Kegg_bar_up <- barplot(kk_up, showCategory=15)+
      labs(title = "Enriched KEGG pathways of up-regulated genes (treat/control)")
    ggsave("results/Kegg_bar_up.pdf", Kegg_bar_up, height = 6, width = 6)


    #################################Down GO KEGG

    genelist_down <- filter(DEGresults, regulation == "Down_reg")$logFC
    names(genelist_down) <- filter(DEGresults, regulation == "Down_reg")$SYMBOL

    GeneList_down <- bitr(names(genelist_down), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    head(GeneList_down)

    gene_down <- GeneList_down$ENTREZID

    ego_all_down <- enrichGO(gene = gene_down,
                             OrgDb = "org.Hs.eg.db",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.5,
                             readable = TRUE)

    ego_bar_down <- barplot(ego_all_down,showCategory=20, split= "ONTOLOGY")+
      facet_grid(ONTOLOGY~.,scale="free")
    # scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    # theme_bw()

    ggsave("results/ego_bar_down.pdf", ego_bar_down, height = 10, width = 10)
    GO_all_down <- as.data.frame(ego_all_down)

    write.csv(GO_all_down, file = "results/GO_all_down.csv")


    kk_down <- enrichKEGG(gene = gene_down,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)
    kk_down <- setReadable(kk_down, keyType = "ENTREZID", OrgDb = "org.Hs.eg.db")

    KEGG_down <- as.data.frame(kk_down)

    write.csv(KEGG_down, file = "results/KEGG_down.csv")

    Kegg_bar_down <- barplot(kk_down, showCategory=15)+
      labs(title = "Enriched KEGG pathways of down-regulated genes (treat/control)")
    ggsave("results/Kegg_bar_down.pdf", Kegg_bar_down, height = 6, width = 6)

    print(list(volcano, top30_up_down, ego_bar_up, ego_bar_down,
               Kegg_bar_up, Kegg_bar_down))
  }
}
