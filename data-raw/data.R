

require(magrittr)
url <- "http://zhongguose.com/colors.json"

ch_color <-  httr::GET(url,encode="json") %>%
  httr::content(.,as = "parsed") %>%
  do.call(rbind, .)%>%
  as.data.frame() %>%
  dplyr::select(-CMYK,-RGB)

save(ch_color, file = "R/ch_color.rda")
