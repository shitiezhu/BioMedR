


qPCR <- function(data = data, ctrl = ctrl, samples = samples, facet = TRUE, filename = filename){

  suppressMessages(library(tidyverse), library(ggpubr))
  data_a <- data %>%
    full_join({group_by(.,Sample.Name) %>%
        summarise_at("Ctrl.CT", ~ mean(.x, na.rm = T))}, by = "Sample.Name", suffix=c("",".mean")) %>%
    mutate(d_Ct = Ct - Ctrl.Ct.mean) %>%
    full_join({filter(.,Sample.Name == ctrl) %>% group_by(Target.Name) %>%
        summarise_at("d_Ct", ~ mean(.x, na.rm = T))}, by = "Target.Name", suffix=c("",".mean")) %>%
    mutate(dd_Ct = d_Ct - d_Ct.mean) %>%
    mutate(`2^-dd_Ct` = 2^-(dd_Ct)) %>%
    mutate(Sample.Name = ordered(Sample.Name, levels = samples))



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
  write.csv(data_a, file.path("results",  paste0(filename,"-qPCR data.csv")), row.names = FALSE)

  ggsave(paste0(filename,"-barplot.jpg"), path = "results", p, width = 8,height = 6)
}

# Sys.time() %>%  stringr::str_extract_all('[0-9]+') %>% .[[1]] %>% paste(collapse = "")



