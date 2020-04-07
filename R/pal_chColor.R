pal_chColor <- function(palette = palette, alpha = 1) {

  palette = sample(as.character(ch_color$hex), 525)
  #palette = as.character(ch_color$hex)
  #palette <- match.arg(palette)

  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1)")

  raw_cols <-  palette
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )

  scales::manual_pal(unname(alpha_cols))
}

scale_color_chColor <- function(palette = palette, alpha = 1, ...) {
  #palette <- match.arg(palette)
  discrete_scale("colour", "chColor", pal_chColor(palette, alpha), ...)
}


scale_colour_chColor <- scale_color_chColor


scale_fill_chColor <- function(palette = palette, alpha = 1, ...) {
  #palette <- match.arg(palette)
  discrete_scale("fill", "chColor", pal_chColor(palette, alpha), ...)
}
