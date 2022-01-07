cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
               "#CC79A7", "#999999", "#FF3300", "#FF99CC", "#3300CC", "#CCFFCC")

#' Plot PCA
plotPCA_pquantr <- function(pdat, meta, with.legend=TRUE, point.size=1.5, text.size=10,
                    label.size=3, legend.text.size=7, palette=cbPalette) {
  pca <- prcomp(t(na.omit(log10(pdat))), scale.=TRUE, center=TRUE)
  
  p <- data.frame(
    x = pca$x[, 1],
    y = pca$x[, 2],
    condition = factor(unique(meta$condition)),
    sample = meta$sample
  )
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  
  g <- ggplot(p, aes_(x=~x, y=~y, label=~sample)) +
    simple_theme +
    theme(
      legend.title = element_text(size = legend.text.size),
      legend.text = element_text(size = legend.text.size)
    ) +
    coord_cartesian(expand=TRUE) +
    geom_point(aes_(colour=~condition), size=point.size) +
    ggrepel::geom_text_repel(size=label.size) +
    scale_color_manual(values=palette) +
    theme(text = element_text(size=text.size)) +
    labs(x=pca1, y=pca2)
  if(!with.legend) g <- g + theme(legend.position="none")
  g
}