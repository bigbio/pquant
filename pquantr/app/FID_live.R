dynamic_fid_selectProtein <- function(tab, input, max_hover=1) {
  sel <- NULL
  tab_idx <- as.numeric(input$allProteinTable_rows_selected)
  if(!is.null(input$plot_brush)){
    #brushed <- na.omit(brushedPoints(tab, input$plot_brush))
    brushed <- brushedPoints(tab, input$plot_brush)
    sel <- as.numeric(rownames(brushed))
  } else if(!is.null(input$plot_hover)) {
    near <- nearPoints(tab, input$plot_hover, threshold = 20, maxpoints = max_hover)
    sel <- as.numeric(rownames(near))
  } else if(length(tab_idx) > 0) {
    sel <- tab_idx
  }
  return(sel)
}

dynamic_fid_replicateTable <- function(tab, input, pdat, max_points) {
  renderTable({
    sel <- dynamic_fid_selectProtein(tab, input)
    if(!is.null(sel)) {
      dat <- pdat[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      if(length(sel) <= max_points) {
        df <- data.frame(
          Sample=colnames(dat),
          Intensity=signif(colMeans(dat, na.rm = TRUE),
                           3))
        df$Intensity[is.nan(df$Intensity)] <- NA
        df
      }
    }
  }, width = "300px")
}

dynamic_fid_significanceTable <- function(tab, res, input, flag) {
  renderDT({
    sel <- dynamic_fid_selectProtein(tab, input)
    if(!is.null(sel)) {
      if(flag == 'n') {
        df <- data.frame(
          `Protein`=sprintf("%s", res$Protein[sel]),
          `GeneID`=sprintf("%s", res$ENTREZID[sel]),
          `GeneName`=sprintf("%s", res$GENENAME[sel]),
          `Log2FC`=sprintf("%.2g", res$log2FC[sel]),
          `P-value`=sprintf("%.2g", res$pvalue[sel]),
          `adjusted P-value`=sprintf("%.2g", res$adj.pvalue[sel]),
          check.names = FALSE)
      } else {
        df <- data.frame(
          `Gene`=sprintf("%s", res$Protein[sel]),
          `GeneName`=sprintf("%s", res$GENENAME[sel]),
          `Log2FC`=sprintf("%.2g", res$log2FC[sel]),
          `P-value`=sprintf("%.2g", res$pvalue[sel]),
          `adjusted P-value`=sprintf("%.2g", res$adj.pvalue[sel]),
          check.names = FALSE)
      }
    return(df)  
    }else { return(NULL) }
  }, width = "100%",
  options = list(lengthChange = FALSE, scrollX = TRUE, searching = FALSE))
}


dynamic_fid_jitterPlot <- function(tab, input, pdat, max_points, meta) {
  renderPlot({
    sel <- dynamic_fid_selectProtein(tab, input)
    if(!is.null(sel) && length(sel) <= max_points) {
      dat <- pdat[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      m <- colMeans(dat, na.rm = TRUE)
      s <- apply(dat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(na.omit(length(x))))
      n <- length(sel)
      p <- data.frame(
        intensity = m,
        lo = m - s,
        up = m + s,
        condition = factor(meta$condition)
      )
      p$shape <- rep(21, length(p$intensity))
      p$shape[which(p$intensity==0)] <- 24
      pd <- ggplot2::position_dodge(width = 0.4)
      ggplot(p, aes_(x=~condition, y=~intensity, ymin=~lo, ymax=~up, shape=~shape)) +
        theme(text = element_text(size=20), legend.position = "none", legend.direction = "horizontal") +
        {if (input$intensityScale == '') ylim(0, NA)} +
        geom_point(position=pd, size=4, na.rm=TRUE) +
        {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
        scale_shape_identity() +  # necessary for shape mapping
        viridis::scale_fill_viridis(discrete=TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
        {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
    }
  })
  
}




simple_theme <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black")
  )
simple_theme_grid <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(colour = "grey95"),
    axis.line = ggplot2::element_line(colour = "black")
  )


dynamic_fid_plotFID <- function(pdat, conditions, pair, bins=80, marginal.histograms=FALSE,
                            xmin=NULL, xmax=NULL, ymax=NULL, text.size=12, point.size=1.6,
                            show.legend=TRUE, plot.grid=TRUE, binhex=TRUE, transform.fun=log10) {
  if(binhex & marginal.histograms) {
    warning("Cannot plot with both binhex=TRUE and marginal.histograms=TRUE. Ignoring binhex.")
    binhex=FALSE
  }
  
  title <- paste(pair, collapse=":")
  
  ttab <- transform.fun(pdat)
  
  # helper function
  
  condMeans <- function(cond) {
    m <- rowMeans(ttab[,which(conditions == cond), drop=FALSE], na.rm=TRUE)
    m[is.nan(m)] <- NA
    m
  }
  
  # build data frame with x-y cooordinates
  # including infinities
  m1 <- condMeans(pair[1])
  m2 <- condMeans(pair[2])
  good <- !is.na(m1) & !is.na(m2)
  df <- data.frame(
    id = rownames(ttab),
    x = (m1 + m2) / 2,
    y = m2 - m1,
    good = good
  )
  
  mx <- 1.1 * max(abs(df$y), na.rm=TRUE)
  m <- rbind(m1[!good], m2[!good])
  df[!good, "x"] <- colSums(m, na.rm=TRUE)
  df[!good, "y"] <- ifelse(is.na(m[1,]), mx, -mx)
  
  g <- ggplot(df[good, ], aes_(x=~x, y=~y))
  
  if(binhex) {
    g <- g + stat_binhex(bins=bins, show.legend=show.legend, na.rm=TRUE) +
      viridis::scale_fill_viridis(name="count", na.value=NA)
    #scale_fill_gradientn(colours=c("seagreen","yellow", "red"), name = "count",na.value=NA)
  } else {
    g <- g + geom_point(size=point.size, na.rm=TRUE) +
      geom_point(data=df[!good,], aes_(x=~x, y=~y), colour="orange", size=point.size, na.rm=TRUE)
  }
  
  if(plot.grid) {
    g <- g + simple_theme_grid
  } else {
    g <- g + simple_theme
  }
  
  g <- g + geom_abline(colour='red', slope=0, intercept=0) +
    labs(title=title, x=paste0(pair[1], '+', pair[2]), y=paste0(pair[2], '-', pair[1])) +
    theme(text = element_text(size=text.size))
  
  if(!is.null(xmin) && !is.null(xmax)) g <- g + scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0))
  if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(-ymax, ymax), expand = c(0, 0))
  if(marginal.histograms) g <- ggExtra::ggMarginal(g, size=10, type = "histogram", xparams=list(bins=100), yparams=list(bins=50))
  return(g)
}