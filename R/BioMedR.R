
add_flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  require(grid)

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels,
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)),
                    length.out = sum(d.select)),
                "npc"))
  }


  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x - unit(0.01, "npc"),
                           x1 = new.label$x + unit(0.1, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4,
                                     l = 4)

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

# heat refers to the original heatmap produced from the pheatmap() function
# kept.labels should be a vector of labels you wish to show
# repel.degree is a number in the range [0, 1], controlling how much the
# labels are spread out from one another


Col2Grey <- function(col){
  rgb <- col2rgb(col)
  g <- rbind( c(0.1, 0.3, 0.6) ) %*% rgb
  rgb(g, g, g, maxColorValue=255)
}

library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

TCGA_KIRC <- data.table::fread("GDC_TCGA_KIRC_YAPgenes.tsv") %>% as.data.frame() #%>%
# mutate_at(7:ncol(.), ~ 2^.x -1) %>%
# na.omit() %>%
# mutate_at(7:ncol(.), fpkmToTpm) %>%
# mutate_at(7:ncol(.), ~log2(.x +1))

SYMBOL = clusterProfiler::bitr(str_sub(colnames(TCGA_KIRC[,7:22]), 1,15), "ENSEMBL", "SYMBOL", OrgDb = "org.Hs.eg.db")

names(TCGA_KIRC)[7:22] <- SYMBOL$SYMBOL

TCGA_KIRC <- TCGA_KIRC %>% dplyr::select(-GPX6)


OSS <- function(data, gene){

  library(tidyverse)
  library(survival)
  library(survminer)
  library(patchwork)

  survdata <- data %>%
    filter( sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
    na.omit() %>%
    mutate(group = ifelse(.[[gene]] > quantile(.[[gene]], 0.5), "high",
                          ifelse(.[[gene]] < quantile(.[[gene]], 0.5),"low", NA)))
  group = survdata$group
  surv <- survdata %>% ggsurvplot(surv_fit(survival::Surv(OS.time/30, OS) ~ group, data = .), data = ., pval = T,
                                  pval.coord = c(0.25, 0.25), size = 0.5, censor.size= 3) %>%
    {.$plot +
        labs(x = "Time (months)")+
        ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
        scale_color_manual(labels = c(paste0("high ", "(n = ", length(which(group == "high")), ")"),
                                      paste0("low  ", "(n = ", length(which(group == "low")), ")")),
                           values = c("#C00750", "#372F82"))+
        labs(title = paste0(gene, " survival plot"))+
        theme(legend.position = c(0.75, 0.90), plot.title = element_text(size = 12))
    }

  N_T <- data %>%
    mutate(type = ifelse(sample_type.samples == "Solid Tissue Normal", "normal", "tumor")) %>%
    na.omit() %>%
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    ggplot(aes(x = genename, y = expression, fill = type))+
    geom_boxplot()+
    ggpubr::stat_compare_means(aes(x = genename, y = expression, group = type,
    ), method = "t.test")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    scale_fill_manual(values = c("#372F82", "#C00750"))+
    theme(legend.position = "right")

  stage <- data %>%
    filter(sample_type.samples != "Solid Tissue Normal",
           tumor_stage.diagnoses != "not reported",
           tumor_stage.diagnoses != "stage x") %>%
    mutate(Stage = str_replace_all(tumor_stage.diagnoses, c("ia" = "i", "ib" = "i", "ic" = "i"))) %>%
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    na.omit() %>%
    ggplot(aes(x = genename,
               y = expression, fill = Stage))+
    geom_boxplot(outlier.size = 1)+
    #ggpubr::stat_compare_means(aes(group = Stage),method = "kruskal.test", label.y = 26)+
    #
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    theme(legend.position = "right")

  ######################################################

  surv_stage <- data %>%
    filter(sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
    na.omit() %>%
    group_by(tumor_stage.diagnoses) %>%
    nest() %>%
    arrange(tumor_stage.diagnoses) %>%
    mutate(sp = map(data, function(df) {
      categro <- surv_cutpoint(time = "OS.time", event = "OS",variables = gene, data = df) %>%
        surv_categorize(labels = c("low", "high"))
      group = categro[,gene]
      surv_plot <- ggsurvplot(surv_fit(survival::Surv(OS.time/30, OS) ~ group, data = categro),
                              data = categro, pval = T, pval.coord = c(0.25, 0.25),
                              size = 0.5, censor.size= 3) %>%
        {.$plot +
            labs(x = "Time (months)")+
            ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
            scale_color_manual(labels = c(paste0("high ", "(n = ", length(which(group == "high")), ")"),
                                          paste0("low  ", "(n = ", length(which(group == "low")), ")")),
                               values = c("#C00750", "#372F82"))+
            theme(legend.position = "bottom", plot.title = element_text(size = 12))}
    }
    )) %>%
    pull(sp) %>%
    wrap_plots(nrow = 1)

  ##############################################

  p <- (surv + N_T + stage)/ surv_stage

  if(!dir.exists("batchplots")) {dir.create("batchplots")}

  ggsave(paste0(gene," plot.jpg"), p, path = "batchplots",dpi = 600, width = 14, height = 8)

  # export::graph2ppt(p, paste0("batchplots/",gene," plot in TCGA_KIRC.pptx"), width = 14, height = 8)

}




