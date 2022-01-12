
Get2019IF <- function(){
  library(httr)
  library(magrittr)

  url <- "https://jcr.clarivate.com/JournalHomeGridJson.action?_dc=1593424942122&jcrYear=2019&edition=Both&categoryIds=&subjectCategoryScheme=WoS&jifQuartile=&impactFactorRangeFrom=&impactFactorRangeTo=&averageJifPercentileRangeFrom=&averageJifPercentileRangeTo=&OAFlag=N&start=0&limit=12855&sort=%5B%7B%22property%22%3A%22journalImpactFactor%22%2C%22direction%22%3A%22DESC%22%7D%5D"
  POST(url, encode = "json") %>%

    content("text", encoding = "UTF-8") %>%
    rjson::fromJSON() %>%
    .$data %>%
    sapply(.,c) %>%
    t() %>%
    as.data.frame()
}

data <- Get2019IF()

openxlsx::write.xlsx(data, file = "C:\\Users\\as\\Desktop\\2019IF.xlsx" )
