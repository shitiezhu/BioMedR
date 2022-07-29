qPCR <- function(data = data, Samples = Samples, facet = TRUE){

  suppressMessages(library(tidyverse), library(ggpubr))
  data_a <- data %>%
    full_join({group_by(.,Sample.Name) %>%
        summarise_at("Ctrl.CT", ~ mean(.x, na.rm = T))}, by = "Sample.Name", suffix=c("",".mean")) %>%
    mutate(d_CT = CT - Ctrl.CT.mean) %>%
    full_join({filter(.,Sample.Name == Ctrl) %>% group_by(Target.Name) %>%
        summarise_at("d_CT", ~ mean(.x, na.rm = T))}, by = "Target.Name", suffix=c("",".mean")) %>%
    mutate(dd_CT = d_CT - d_CT.mean) %>%
    mutate(`2^-dd_CT` = 2^-(dd_CT)) %>%
    mutate(Sample.Name = ordered(Sample.Name, levels = Samples))



  # install.packages("ggpubr")
  p <- ggbarplot(data_a, x = "Target.Name", y = "2^-dd_CT", add = "mean_se",
            color = "Sample.Name", palette = "jco", fill = "Sample.Name",
            position = position_dodge(0.8))+
    stat_compare_means(aes(group = Sample.Name), label = "p.format") +
    labs(y = "Rel. mRNA expression", x = "")
  if (facet == TRUE){
   p <- p + facet_wrap(~Target.Name, scales = "free")
  }


  if(!dir.exists("results")) {dir.create("results")}
  write.csv(data_a, "results/qPCR data ready for analysis.csv", row.names = FALSE)

  ggsave("results/barplot.jpg", p, width = 8,height = 6)
}
