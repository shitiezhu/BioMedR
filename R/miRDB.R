miRDB <- function(searchBox = NA, searchType = c("miRNA", "gene"), Species = c("Human", "Mouse")){

  suppressWarnings(suppressPackageStartupMessages({
    library(httr)
    library(curl)
    library(rlist)
    library(tidyverse)
    library(rvest)} ))

  headers <- c('Accept'='text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
               'Content-Type'='application/x-www-form-urlencoded',
               'User-Agent'='Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.88 Safari/537.36',
               'Referer'='http://mirdb.org/index.html',
               'Connection'='keep-alive'
  )
  payload<-list(
    Species = match.arg(Species),
    searchBox = searchBox,
    searchType = match.arg(searchType),
    submitButton = "Go",
    geneChoice = "geneID"
  )

  url <- "http://mirdb.org/cgi-bin/search.cgi"

  POST(url,add_headers(.headers =headers),body = payload, encode= "form" ) %>%
    read_html() %>%
    html_nodes("table#table1 tr") %>% html_text() %>% str_split("\\n") %>%
    do.call(rbind,.) %>% .[-1,] %>%
    data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("Target Detail", "Target Rank", "Target Score", "miRNA Name", "Gene Symbol", "Gene Description")) %>%
    .[,-7]
}

df <- miRDB(searchBox = "hsa-miR-145-3p",
            searchType = "miRNA",
            Species = "Human")

df2 <- miRDB(searchBox = "7157",
             searchType = "gene",
             Species = "Human")
