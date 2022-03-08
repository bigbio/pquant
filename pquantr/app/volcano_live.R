dynamic_selectProtein <- function(tab, input, max_hover=1) {
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

dynamic_replicateTable <- function(tab, input, pdat, max_points) {
  renderTable({
    sel <- dynamic_selectProtein(tab, input)
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

dynamic_significanceTable <- function(tab, res, input, flag) {
  renderDT({
    sel <- dynamic_selectProtein(tab, input)
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
    } else { return(NULL) }
  }, width = "100%",
  options = list(lengthChange = FALSE, scrollX = TRUE, searching = FALSE))
}


dynamic_jitterPlot <- function(tab, input, pdat, max_points, meta) {
  renderPlot({
    sel <- dynamic_selectProtein(tab, input)
    if(!is.null(sel) && length(sel) <= max_points) {
      dat <- pdat[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      m <- colMeans(dat, na.rm = TRUE)
      s <- apply(dat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(na.omit(length(x))))
      n <- length(sel)
      
      meta_conditions <- meta$condition
      for (i in 1:length(meta_conditions)){
        
          meta_cond <- meta_conditions[i]
          
          if (nchar(meta_cond) > 32) {
              cond <- unlist(strsplit(meta_cond, split = " "))
              
              tmp_abb = NULL
          
              for (str in cond) {
                abb <- substr(str, 0, 1)  #name is too long, only to extract the first letter of each word
                tmp_abb <- paste(tmp_abb, abb, sep = " ")
              }
              
              meta_conditions[i] <- substr(tmp_abb, 2, nchar(tmp_abb))
          }
      }
      
      p <- data.frame(
        intensity = m,
        lo = m - s,
        up = m + s,
        condition = factor(meta_conditions)
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



dynamic_plotVolcano <- function(res, bins=80, xmax=NULL, ymax=NULL, marginal.histograms=FALSE, text.size=12, show.legend=TRUE,
                                plot.grid=TRUE, binhex=TRUE) {
  if(binhex & marginal.histograms) {
    warning("Cannot plot with both binhex=TRUE and marginal.histograms=TRUE. Ignoring binhex.")
    binhex=FALSE
  }
  
  tr <- attr(res, "transform.fun")
  conds <- attr(res, "conditions")
  xlab <- ifelse(is.null(tr), "log2FC", paste(tr, "log2FC"))
  tit <- paste(conds, collapse=":")
  id <- names(res)[1]
  
  g <- ggplot(res, aes_(~log2FC, ~-log10(pvalue)))
  
  if(binhex) {
    g <- g + stat_binhex(bins=bins, show.legend=show.legend, na.rm=TRUE) +
      viridis::scale_fill_viridis(name="count", na.value=NA)
    #scale_fill_gradientn(colours=c("seagreen","yellow", "red"), name = "count", na.value=NA)
  } else {
    g <- g + geom_point(na.rm=TRUE)
  }
  
  if(plot.grid) {
    g <- g + simple_theme_grid
  } else {
    g <- g + simple_theme
  }
  
  g <- g + geom_vline(colour='red', xintercept=0) +
    theme(text = element_text(size=text.size)) +
    labs(x=xlab, y="-log10 p", title=tit)
  
  
  if(!is.null(xmax)) g <- g + scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0))
  if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0))
  
  if(marginal.histograms) g <- ggExtra::ggMarginal(g, size=10, type = "histogram", xparams=list(bins=100), yparams=list(bins=50))
  return(g)
}
