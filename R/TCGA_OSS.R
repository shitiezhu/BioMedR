OSS <- function(data = data, gene = gene, cut.point = "0.5", out.format = "r"){

  library(tidyverse)
  library(survival)
  library(survminer)
  library(patchwork)
  data = data
  gene = gene

  if(cut.point == "0.5")  {
    survdata <- data %>%
      filter( sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
      mutate(group = ifelse(.[[gene]] > quantile(.[[gene]], 0.5), "high",
                            ifelse(.[[gene]] < quantile(.[[gene]], 0.5),"low", NA)))
  } else if(cut.point == "0.75"){
    survdata <- data %>%
      filter( sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
      mutate(group = ifelse(.[[gene]] > quantile(.[[gene]], 0.75), "high",
                            ifelse(.[[gene]] < quantile(.[[gene]], 0.75),"low", NA)))
  } else if(cut.point == "0.75/0.25"){
    survdata <- data %>%
      filter( sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
      mutate(group = ifelse(.[[gene]] > quantile(.[[gene]], 0.75), "high",
                            ifelse(.[[gene]] < quantile(.[[gene]], 0.25),"low", NA)))
  } else if(cut.point == "best"){
    survdata <- data %>%
      filter( sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
      mutate(group = {surv_cutpoint(time = "OS.time", event = "OS",variables = gene, data = .) %>%
          surv_categorize(labels = c("low", "high")) %>% pull(gene)})
  }

  group = survdata$group
  surv <- survdata %>% ggsurvplot(surv_fit(survival::Surv(OS.time/30, OS) ~ group, data = ., check.names = T), data = ., pval = T,
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
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    ggplot(aes(x = genename, y = expression, color = type))+
    geom_boxplot()+
    geom_point(position = position_jitterdodge(0.15),size  = 1) +
    ggpubr::stat_compare_means(aes(x = genename, y = expression, group = type,
    ), method = "t.test")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    scale_color_manual(values = c("#372F82", "#C00750"))+
    theme(legend.position = "right")

  ######################################################
  stage1 <- data %>%
    filter(sample_type.samples != "Solid Tissue Normal",
           tumor_stage.diagnoses != "not reported",
           tumor_stage.diagnoses != "stage x") %>%
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    ggplot(aes(x = genename,
               y = expression, color = tumor_stage.diagnoses))+
    geom_boxplot(outlier.size = 1)+
    geom_point(position = position_jitterdodge(0.1),size  = 1) +
    #ggpubr::stat_compare_means(aes(group = Stage),method = "kruskal.test", label.y = 26)+
    #
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    theme(legend.position = "right")


  surv_stage1 <- data %>%
    filter(sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
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

  ##############################

  stage2 <- data %>%
    filter(sample_type.samples != "Solid Tissue Normal",
           tumor_stage.diagnoses != "not reported",
           tumor_stage.diagnoses != "stage x") %>%
    mutate(stage = ifelse(tumor_stage.diagnoses == "stage i" | tumor_stage.diagnoses == "stage ii",
                          "stage I + II", "stage III + IV")) %>%
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    ggplot(aes(x = genename,
               y = expression, color = stage))+
    geom_boxplot(outlier.size = 1)+
    geom_point(position = position_jitterdodge(0.1),size  = 1) +
    ggpubr::stat_compare_means(aes(x = genename, y = expression, group = stage), method = "t.test")+
    #ggpubr::stat_compare_means(aes(group = Stage),method = "kruskal.test", label.y = 26)+
    #
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    theme(legend.position = "right")

  surv_stage2 <- data %>%
    filter(sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
    mutate(stage = ifelse(tumor_stage.diagnoses == "stage i" | tumor_stage.diagnoses == "stage ii", "stage I + II", "stage III + IV")) %>%
    group_by(stage) %>%
    nest() %>%
    arrange(stage) %>%
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
  ################################
  stage3 <- data %>%
    filter(sample_type.samples != "Solid Tissue Normal",
           tumor_stage.diagnoses != "not reported",
           tumor_stage.diagnoses != "stage x") %>%
    mutate(stage = ifelse(tumor_stage.diagnoses == "stage i","stage I", "stage II + III + IV")) %>%
    pivot_longer(gene, names_to = "genename", values_to = "expression") %>%
    ggplot(aes(x = genename,
               y = expression, color = stage))+
    geom_boxplot(outlier.size = 1)+
    geom_point(position = position_jitterdodge(0.1),size  = 1) +
    ggpubr::stat_compare_means(aes(x = genename, y = expression, group = stage), method = "t.test")+
    #ggpubr::stat_compare_means(aes(group = Stage),method = "kruskal.test", label.y = 26)+
    #
    labs(x = "", y = "mRNA expression (log2FPKM)")+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 12)+
    theme(legend.position = "right")

  surv_stage3 <- data %>%
    filter(sample_type.samples == "Primary Tumor", tumor_stage.diagnoses != "not reported") %>%
    mutate(stage = ifelse(tumor_stage.diagnoses == "stage i", "stage I", "stage II + III + IV")) %>%
    group_by(stage) %>%
    nest() %>%
    arrange(stage) %>%
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

  p <- (surv + N_T + stage1)/ surv_stage1 / (surv_stage2 + stage2)/ (surv_stage3 + stage3)

  if(!dir.exists("batchplots")) {dir.create("batchplots")}

  if(out.format == "jpg") {ggsave(paste0(gene," plot in TCGA.jpg"), p, path = "batchplots",dpi = 600, width = 14, height = 16)
  } else if (out.format == "pptx") {
    export::graph2ppt(p, paste0("batchplots/",gene," plot in TCGA.pptx"), width = 14, height = 16)
  }else if (out.format == "pdf") {
    ggsave(paste0(gene," plot in TCGA.pdf"), p, path = "batchplots",dpi = 600, width = 14, height = 16)
  }else if (out.format == "r") {
    p
  }
  #
}
