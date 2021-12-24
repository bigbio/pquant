cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
               "#CC79A7", "#999999", "#FF3300", "#FF99CC", "#3300CC", "#CCFFCC")

#' Plot peptide or protein count per sample
plotCount <- function(pdat, meta, x.text.size=10, palette=cbPalette){
    entry.count <- apply(pdat, 2, function(x) sum(!is.na(x)))
    med.count <- median(entry.count)
    condition = meta$condition
    df <- data.frame(x=meta$sample, y=entry.count, condition=condition)
    g <- ggplot(df, aes_(x=~x,y=~y,fill=~condition)) +
      geom_col(colour='grey60') +
      simple_theme +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=x.text.size)) +
      labs(x="Sample", y="Count") +
      labs(title = paste0("Median count = ", med.count)) +
      theme(plot.title=element_text(hjust=0, size=12))
    if(nlevels(as.factor(df$condition)) <= length(palette)) g <- g + scale_fill_manual(values=palette)
    #g <- g + viridis::scale_fill_viridis(discrete=TRUE)
    g
}


#' Jaccard similarity
jaccardSimilarity <- function(x, y) {
  stopifnot(length(x) == length(y))
  
  intersection <- which(!is.na(x) & !is.na(y))
  union <- which(!is.na(x) | !is.na(y))
  jaccard <- ifelse(length(union) > 0, length(intersection) / length(union), 0)
  return(jaccard)
}

#' Detection Jaccard similarity
plotDetectionSimilarity <- function(pdat, bin.size=0.01, text.size=12, plot.grid=TRUE, hist.colour='grey30') {
  
  n <- ncol(pdat)
  # indices of upper triangular: all pair-wise combinations of columns in tab
  pairs <- which(upper.tri(matrix(1, n, n)) == TRUE, arr.ind=TRUE)
  # Jaccard simliarity for each pair of columns in tab
  sim <- apply(pairs, 1, function(i) jaccardSimilarity(pdat[,i[1]], pdat[,i[2]]))
  
  ggplot(data.frame(sim), aes_(~sim, ~..density..)) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    geom_histogram(breaks=seq(0, 1, bin.size), colour=hist.colour, fill=hist.colour) +
    labs(x='Jaccard similarity', y='Density') +
    theme(text = element_text(size=text.size))
}


#' Plot distance matrix
plotDistanceMatrix <- function(pdat, meta, distance=c("correlation"), text.size=10) {
  distance <- match.arg(distance)
  
  corr.mat <- cor(pdat, use="complete.obs")
  m <- reshape2::melt(corr.mat, varnames=c("Sample1", "Sample2"))
  m$Sample1 <- factor(m$Sample1, levels=meta$sample)
  m$Sample2 <- factor(m$Sample2, levels=meta$sample)
  ggplot(m, aes_(x=~Sample1, y=~Sample2)) +
    geom_tile(aes_(fill=~value)) +
    viridis::scale_fill_viridis() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=text.size),
      axis.text.y = element_text(size=text.size)
    ) +
    labs(x='Sample', y='Sample', fill="Correlation")
}


#' Plot clustering dendrogram
plotClustering <- function(pdat, x.text.size=10) {
  corr.mat <- cor(pdat, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  hc <- hclust(dis)
  dendr <- ggdendro::dendro_data(hc)
  dat <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5, size=x.text.size),
    axis.line.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    geom_segment(data=dat, aes_(x=~x, y=~y, xend=~xend, yend=~yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(dat$y) * 1.03)) +
    labs(x="Sample", y="Distance")
}


#' Plot PCA
plotPCA <- function(pdat, meta, with.legend=TRUE, point.size=1.5, text.size=10,
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


#' Plot distribution of intensities/ratios for each sample
plotSampleDistributions <-
  function(pdat, meta, title="", method=c("violin", "dist", "box"), x.text.size=7, n.grid.rows=3,
           hist.bins=100, x.text.angle=90, vmin=as.numeric(NA), vmax=as.numeric(NA), log.scale=TRUE,
           log.base=10, palette=cbPalette, fill=NULL, colour=NULL, hline=FALSE) {
    method <- match.arg(method)
    
    m <- reshape2::melt(pdat, varnames=c("ID", "sample"))
    mt <- data.frame(meta, row.names = meta$sample)
    if(!is.null(fill)) m[['fill']] <- as.factor(mt[m$sample, fill])
    if(!is.null(colour)) m[['colour']] <- as.factor(mt[m$sample, colour])
    
    if(log.scale > 0) m$value <- log(m$value, base=log.base)
    lg <- ifelse(log.scale, paste("log", log.base), "")
    
    if(method == "box" | method == "violin") {
      g <- ggplot(m, aes_(x=~sample, y=~value)) +
        simple_theme +
        ylim(vmin, vmax) +
        theme(axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust=0.5, size=x.text.size)) +
        labs(title=title, x="sample", y=paste0(lg, " value"))
      if(hline) g <- g + geom_hline(yintercept=0, colour='grey')
      if(method=="box") g <- g + geom_boxplot(outlier.shape=NA, na.rm=TRUE)
      if(method=="violin") g <- g + geom_violin(na.rm=TRUE, draw_quantiles=c(0.25, 0.5, 0.75), scale="width")
    } else if (method == "dist") {
      g <- ggplot(m, aes_(x=~value)) +
        geom_histogram(bins=hist.bins) +
        facet_wrap(~sample, nrow=n.grid.rows) +
        xlim(vmin, vmax)
    } else stop("Wrong method.")
    
    if(!is.null(fill)) {
      g <- g + aes_(fill=~fill)
      if(nlevels(m$fill) <= length(palette)) g <- g + scale_fill_manual(name=fill, values=palette)
    }
    if(!is.null(colour)) {
      g <- g + aes_(colour=~colour)
      if(nlevels(m$colour) <= length(palette)) g <- g + scale_color_manual(name=colour, values=palette)
    }
    g
  }


#' Statistics for an intensity table
intensityStats <- function(pdat, meta) {
  stats <- NULL
  for(cond in unique(meta$condition)) {
    w <- pdat[,which(factor(unique(meta$condition)) == cond), drop=FALSE]
    m <- rowMeans(w, na.rm=TRUE)
    m[which(is.nan(m))] <- NA
    v <- apply(w, 1, function(v) sd(v, na.rm=TRUE)^2)
    ngood <- apply(w, 1, function(v) sum(!is.na(v)))
    stats <- rbind(stats, data.frame(id=rownames(w), condition=cond, mean=m, variance=v, ngood=ngood))
    rownames(stats) <- NULL
  }
  return(stats)
}

#' Plot mean-variance relationship
plotMV <- function(pdat, meta, with.loess=FALSE, bins=80, xmin=5, xmax=10, ymin=7, ymax=20, with.n=FALSE, mid.gradient=0.3, text.size=10) {
  stats <- intensityStats(pdat, meta)
  stats <- stats[which(!is.na(stats$mean) & !is.na(stats$variance) & stats$mean > 0 & stats$variance > 0),]
  stats$mean <- log10(stats$mean)
  stats$variance <- log10(stats$variance)
  protnum <- as.data.frame(table(stats$condition))  # number of proteins in each condition
  colnames(protnum) <- c("condition", "n")
  protnum$labn <- paste0("n = ", protnum$n)
  
  # has to be calculated for each condition separately
  if(with.loess) {
    ldf <- NULL
    for(cond in meta$condition)
    {
      st <- stats[which(stats$condition == cond),]
      ls <- loess(variance ~ mean, data=st)
      x <- seq(from=min(na.omit(st$mean)), to=max(na.omit(st$mean)), by=0.01)
      pr <- predict(ls, x)
      ldf <- rbind(ldf, data.frame(condition=cond, x=x, y=pr))
    }
  }
  
  g <- ggplot(stats, aes_(x=~mean, y=~variance)) +
    simple_theme +
    theme(panel.border = element_rect(fill=NA, color='black')) +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    facet_wrap(~condition) +
    stat_binhex(bins=bins) +
    theme(text = element_text(size=text.size)) +
    viridis::scale_fill_viridis(values=c(0, mid.gradient, 1), name="count", na.value=NA) +
    #scale_fill_gradientn(colours=c("seagreen","yellow", "red"), values=c(0, mid.gradient, 1), name="count", na.value=NA)
    if(with.n) g <- g + geom_text(data=protnum, aes_(x=~xmin+0.5, y=~ymax, label=~labn))
  if(with.loess) g <- g + geom_line(data=ldf, aes_(x=~x,y=~y), color='black')
  return(g)
}