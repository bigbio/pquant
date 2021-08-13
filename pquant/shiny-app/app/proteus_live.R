proteus_selectProtein <- function(tab, input, max_hover=1) {
  sel <- NULL
  tab_idx <- as.numeric(input$allProteinTable_rows_selected)
  if(!is.null(input$plot_brush)){
    brushed <- na.omit(brushedPoints(tab, input$plot_brush))
    sel <- as.numeric(rownames(brushed))
  } else if(!is.null(input$plot_hover)) {
    near <- nearPoints(tab, input$plot_hover, threshold = 20, maxpoints = max_hover)
    sel <- as.numeric(rownames(near))
  } else if(length(tab_idx) > 0) {
    sel <- tab_idx
  }
  return(sel)
}

proteus_replicateTable <- function(tab, input, pdat, max_points) {
  renderTable({
    sel <- proteus_selectProtein(tab, input)
    if(!is.null(sel)) {
      dat <- pdat$tab[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      if(length(sel) <= max_points) {
        df <- data.frame(Sample=colnames(dat), Intensity=signif(colMeans(dat, na.rm = TRUE), 3))
        df$Intensity[is.nan(df$Intensity)] <- NA
        df
      }
    }
  }, width = "80px")
}

proteus_significanceTable <- function(tab, res, input) {
  renderTable({
    sel <- proteus_selectProtein(tab, input)
    if(!is.null(sel) && length(sel) == 1) {
      data.frame(`P-value`=sprintf("%.2g", res$P.Value[sel]), `adjusted P-value`=sprintf("%.2g", res$adj.P.Val[sel]), check.names = FALSE)
    }
  }, width = "100px")
}


proteus_jitterPlot <- function(tab, input, pdat, max_points) {
  renderPlot({
    sel <- proteus_selectProtein(tab, input)
    if(!is.null(sel) && length(sel) <= max_points) {
      dat <- pdat$tab[sel,, drop=FALSE]
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
        condition = factor(pdat$metadata$condition, levels=pdat$conditions),
        replicate = as.factor(pdat$metadata$replicate)
      )
      p$shape <- rep(21, length(p$intensity))
      p$shape[which(p$intensity==0)] <- 24
      pd <- ggplot2::position_dodge(width = 0.4)
      ggplot(p, aes_(x=~condition, y=~intensity, ymin=~lo, ymax=~up, colour=~replicate, shape=~shape, fill=~replicate)) +
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



proteus_allProteinTable <- function(res) {
  DT::renderDataTable({
    # assume first column is id ("protein" or "peptide")
    idcol <- names(res)[1]
    cols <- c(idcol, "logFC", "adj.P.Val", grep("mean_", colnames(res), value=TRUE))
    d <- res[, cols]
    d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) signif(x, 3))
    d <- DT::datatable(d, class = 'cell-border strip hover')
    DT::formatStyle(d, 0, cursor = 'pointer')
  })
}




proteus_proteinInfo <- function(tab, input, pdat, max_points) {
  renderTable({
    sel <- proteus_selectProtein(tab, input)
    if(!is.null(sel)) {
      n <- length(sel)
      if (is.null(pdat$annotation)){
        data.frame(Error='No annotation found on the Proteus object. Consult vignette.')
      } else {
        if (n >= 1 && n <= max_points && sel > 0) {
          data.frame(pdat$annotation[sel, ])
        } else if (n > max_points) {
          data.frame(Error=paste('Only', max_points, 'points can be selected.'))
        }
      }
    }
  })
}