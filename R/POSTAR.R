

POSTAR <- function(key = NA, method = "eclip", protocol = ''){

  library(httr)
  library(jsonlite)
  library(curl)
  library(rlist)
  library(UpSetR)
  library(tidyverse)
  library(rjson)

  headers <- c('Accept'='application/json, text/javascript, */*; q=0.01',
               'Content-Type'='application/x-www-form-urlencoded; charset=UTF-8',
               'User-Agent'='Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36',
               'Referer'='http://lulab.life.tsinghua.edu.cn/postar/gene2.php',
               'Connection'='keep-alive'
  )
  payload<-list(
    key= key,
    searchID= 'searchRBPname',
    method= method,
    species= 'human',
    chr= 'chr1',
    start= '',
    end= '',
    inputState= 2,
    protocol= protocol
  )

  url <- "http://lulab.life.tsinghua.edu.cn/postar/script/ajax/gene.php"

  POST(url,add_headers(.headers =headers),body = payload, encode= "form" ) %>%
    httr::content(., as = "text") %>%
    rjson::fromJSON() %>%
    pluck("data") %>%
    as.data.frame() %>%
    setNames(key)
  #%>%
  # separate_rows(.,.,sep = ",") %>%
  # apply( 1, function(x){
  #   str_extract(x, '[0-9A-Z]+')}
  # ) %>%
  # .[-1] %>%
  # head(-1)%>%
  # as.data.frame(stringsAsFactors = FALSE) %>%
  # plyr::rename(c("." = key))
  # library("rjson")
}
