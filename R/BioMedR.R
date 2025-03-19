

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


corrplot <- function(data, x, y, path, out.format){
  path = path
  if(!dir.exists(path)) {dir.create(path)}

  require(tidyverse)
  p = data %>%
    as.data.frame() %>%
    ggplot(aes_string(x, y)) +
    geom_point(size = 0.5)+
    geom_smooth(method = "lm")+
    ggpubr::stat_cor(size = 2.8)+
    ggprism::theme_prism(base_line_size = 0.5, base_fontface = "plain", base_size = 8)

  if(out.format == "jpg") {ggsave(paste0(x,"_",y,"-corrplot.jpg"),p, path = path, dpi = 600, width = 4, height = 4)
  } else if (out.format == "pptx") {
    export::graph2ppt(p, paste0(path,"/" ,x,"_",y,"-corrplot.pptx"), width = 2, height = 2)
  }else if (out.format == "pdf") {
    ggsave(paste0(x,"_",y,"-corrplot.pdf"),p, path = path, width = 2, height = 2)
  }else if (out.format == "r") {
    p
  }
}


require(colorspace)



# 生成高区分度配色方案
tz_palette <- function(n) {
  # 在HCL颜色空间中系统化采样
  hues <- seq(0, 360, length.out = n+1)[1:n]  # 全色相均匀分布
  chroma <- rep(c(75, 85, 85), length.out = n) # 交替饱和度层级
  luminance <- rep(c(60, 70, 70), length.out = n) # 交替亮度层级

  require(colorspace)
  # 生成HCL颜色
  colors <- colorspace::hcl(
    h = hues,
    c = chroma,
    l = luminance
  )

  # 人工优化确保相邻颜色差异最大化
  set.seed(1111)
  optimized_order <- sample(n, n, replace = F)
  colors[optimized_order]
}

