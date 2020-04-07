
Col2Grey <- function(col){
  rgb <- col2rgb(col)
  g <- rbind( c(0.1, 0.3, 0.6) ) %*% rgb
  rgb(g, g, g, maxColorValue=255)
}
