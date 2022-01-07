library(ggrepel)

groupComparisonPlots_TMT = function(
  data, type, sig = 0.05, FCcutoff = FALSE, logBase.pvalue = 10, ylimUp = FALSE,
  ylimDown = FALSE, xlimUp = FALSE, x.axis.size = 10, y.axis.size = 10, 
  dot.size = 3, text.size = 4, text.angle = 0, legend.size = 13, 
  ProteinName = TRUE, colorkey = TRUE, numProtein = 100, clustering = "both", 
  width = 10, height = 10, which.Comparison = "all", which.Protein = "all",
  address = ""
) {
  Label = Protein = NULL
  
  type = toupper(type)
  input = data.table::as.data.table(data)
  all_labels = as.character(unique(data$Label))
  log_base_FC = ifelse(is.element("log2FC", colnames(data)), 2, 10)
  
  chosen_labels = .checkGCPlotsInput(type, logBase.pvalue, which.Comparison,
                                     all_labels)
  input = input[Label %in% chosen_labels]
  input[, Protein := factor(Protein)]
  input[, Label := factor(Label)]
  
  if (type == "HEATMAP") {
    .plotHeatmap(input, logBase.pvalue, ylimUp, FCcutoff, sig, clustering, 
                 numProtein, colorkey, width, height, log_base_FC,
                 x.axis.size, y.axis.size, address)
  }
  if (type == "VOLCANOPLOT") {
    .plotVolcano(input, which.Comparison, address, width, height, logBase.pvalue,
                 ylimUp, ylimDown, FCcutoff, sig, xlimUp, ProteinName, dot.size,
                 text.size, legend.size, x.axis.size, y.axis.size, log_base_FC)
  }
  if (type == "COMPARISONPLOT") {
    .plotComparison(input, which.Protein, address, width, height, sig, ylimUp, 
                    ylimDown, text.angle, dot.size, x.axis.size, y.axis.size,
                    log_base_FC)
  }
}

#' Preprocess data for volcano plots and create them
#' @inheritParams groupComparisonPlots
#' @keywords internal
.plotVolcano = function(
  input, which.Comparison, address, width, height, log_base_pval,
  ylimUp, ylimDown, FCcutoff, sig, xlimUp, ProteinName, dot.size,
  text.size, legend.size, x.axis.size, y.axis.size, log_base_FC
) {
  adj.pvalue = colgroup = logFC = Protein = issue = Label = newlogFC = NULL
  
  log_adjp = paste0("log", log_base_pval, "adjp")
  all_labels = unique(input$Label)
  input = input[!is.na(adj.pvalue), ]
  colname_log_fc = intersect(colnames(input), c("log2FC", "log10FC"))
  data.table::setnames(input, colname_log_fc, c("logFC"))
  
  if (address == FALSE) {
    if (which.Comparison == 'all') {
      if (length(unique(input$Label)) > 1) {
        stop('** Cannnot generate all volcano plots in a screen. Please set one comparison at a time.')
      }
    } else if (length(which.Comparison) > 1) {
      stop( '** Cannnot generate multiple volcano plots in a screen. Please set one comparison at a time.' )
    }
  }
  
  if (is.numeric(ylimUp)) {
    y.limUp = ylimUp 
  } else {
    y.limUp = ifelse(log_base_pval == 2, 30, 10)
  }
  input[, adj.pvalue := ifelse(adj.pvalue < log_base_pval ^ (-y.limUp),
                               log_base_pval ^ (-y.limUp), adj.pvalue)]
  
  if (!FCcutoff) { 
    logFC_cutoff = 0
  } else {
    logFC_cutoff = log(FCcutoff, log_base_FC)
  }
  input[, colgroup := ifelse(adj.pvalue >= sig, "black",
                             ifelse(logFC > logFC_cutoff,
                                    "red", "blue"))]
  input[, colgroup := factor(colgroup, levels = c("black", "blue", "red"))]
  input[, Protein := as.character(Protein)]
  input[!is.na(issue) & issue == "oneConditionMissing", 
        Protein := paste0("*", Protein)]
  
  savePlot(address, "VolcanoPlot", width, height)
  for (i in seq_along(all_labels)) {
    label_name = all_labels[i]
    single_label = input[Label == label_name, ]
    
    y.limup = ceiling(max(-log(single_label[!is.na(single_label$adj.pvalue), "adj.pvalue"], log_base_pval)))
    if (y.limup < (-log(sig, log_base_pval))) {
      y.limup = (-log(sig, log_base_pval) + 1) ## for too small y.lim
    }
    y.limdown = ifelse(is.numeric(ylimDown), ylimDown, 0)
    x_ceiling = ceiling(max(abs(single_label[!is.na(single_label$logFC) & is.finite(single_label$logFC), logFC])))
    x.lim = ifelse(is.numeric(xlimUp), xlimUp, ifelse((x_ceiling < 3), 3, x_ceiling))
    
    single_label[[log_adjp]] = -log(single_label$adj.pvalue, log_base_pval)
    single_label$newlogFC = single_label$logFC
    single_label[!is.na(issue) &
                   issue == "oneConditionMissing" & 
                   logFC == Inf, newlogFC := (x.lim - 0.2)]
    single_label[!is.na(issue) & 
                   issue == "oneConditionMissing" & 
                   logFC == (-Inf), newlogFC := (x.lim - 0.2) * (-1)]
    plot = .makeVolcano(single_label, label_name, log_base_FC, log_base_pval, x.lim, ProteinName, dot.size,
                        y.limdown, y.limup, text.size, FCcutoff, sig, x.axis.size, y.axis.size,
                        legend.size, log_adjp)
    print(plot)
  }
  if (address != FALSE) {
    dev.off()
  }
}


.checkGCPlotsInput = function(type, log_base, selected_labels, all_labels) {
  checkmate::assertChoice(type, c("HEATMAP", "VOLCANOPLOT", "COMPARISONPLOT"))
  checkmate::assertChoice(log_base, c(2, 10))
  if (selected_labels != "all") {
    if (is.character(selected_labels)) {
      chosen_labels = selected_labels
      wrong_labels = setdiff(chosen_labels, all_labels)
      if (length(wrong_labels) > 0) {
        msg_1 = paste("Please check labels of comparisons.",
                      "Result does not have the following comparisons:")
        msg_2 = paste(wrong_labels, sep = ", ", collapse = ", ")
        msg = paste(msg_1, msg_2)
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
      }
    }
    if (is.numeric(selected_labels)) {
      n_labels = length(all_labels)
      if (n_labels < max(selected_labels)) {
        msg = paste("Please check your selection of comparisons. There are",
                    n_labels, "comparisons in this result.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
      } else {
        chosen_labels = all_labels[selected_labels]
      }
    }  
  } else {
    chosen_labels = all_labels
  }
  chosen_labels
}


.makeVolcano = function(
  input, label_name, log_base_FC, log_base_pval, x.lim, ProteinName, dot.size,
  y.limdown, y.limup, text.size, FCcutoff, sig, x.axis.size, y.axis.size,
  legend.size, log_adjp
) {
  Protein = NULL
  
  plot = ggplot(aes_string(x = "logFC", 
                           y = log_adjp,
                           color = "colgroup",
                           label = "Protein"),
                data = input) +
    geom_point(size = dot.size) +
    scale_colour_manual(values = c("gray65", "blue", "red"), 
                        limits = c("black", "blue", "red"), 
                        breaks = c("black", "blue", "red"), 
                        labels = c("No regulation", "Down-regulated", "Up-regulated")) +
    scale_y_continuous(paste0("-Log", log_base_pval, " (adjusted p-value)"),
                       limits = c(y.limdown, y.limup)) +
    labs(title = unique(label_name))
  plot = plot +
    scale_x_continuous(paste0("Log", log_base_pval, " fold change"),
                       limits = c(-x.lim, x.lim))
  if (ProteinName) {
    if (!(length(unique(input$colgroup)) == 1 & any(unique(input$colgroup) == "black"))) {
      plot = plot +
        geom_text_repel(data = input[input$colgroup != "black", ],
                        aes(label = Protein),
                        size = text.size,
                        col = "black")
    }
  } 
  if (!FCcutoff | is.numeric(FCcutoff)) {
    l = ifelse(!FCcutoff, 20, 10)
    sigcut = data.table::setnames(
      data.table::data.table("sigline", 
                             seq(-x.lim, x.lim, length.out = l),
                             (-log(sig, base = log_base_pval)),
                             "twodash"), 
      c("Protein", "logFC", log_adjp, "line"))
  }
  if (!FCcutoff) {
    plot = plot +
      geom_line(data = sigcut,
                aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                colour = "darkgrey",
                size = 0.6,
                show.legend = TRUE) +
      scale_linetype_manual(values = c("twodash" = 6),
                            labels = c(paste0("Adj p-value cutoff (", sig, ")"))) +
      guides(colour = guide_legend(override.aes = list(linetype = 0)),
             linetype = guide_legend())
  }
  if (is.numeric(FCcutoff)) {
    FCcutpos = data.table::setnames(data.table("sigline", 
                                               log(FCcutoff, log_base_pval), 
                                               seq(y.limdown, y.limup, length.out = 10), 
                                               "dotted"),
                                    c("Protein", "logFC", log_adjp, "line"))
    FCcutneg = data.table::setnames(data.table("sigline", 
                                               (-log(FCcutoff, log_base_pval)), 
                                               seq(y.limdown, y.limup, length.out = 10), 
                                               "dotted"),
                                    c("Protein", "logFC", log_adjp, "line"))
    plot = plot +
      geom_line(data = sigcut, 
                aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                colour = "darkgrey",
                size = 0.6,
                show.legend = TRUE) +
      geom_line(data = FCcutpos,
                aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                colour = "darkgrey",
                size = 0.6,
                show.legend = TRUE) +
      geom_line(data = FCcutneg,
                aes_string(x = "logFC", y = log_adjp, linetype = "line"),
                colour = "darkgrey",
                size = 0.6) +
      scale_linetype_manual(values = c("dotted" = 3, "twodash" = 6),
                            labels = c(paste0("Fold change cutoff (", FCcutoff, ")"),
                                       paste0("Adj p-value cutoff (", sig, ")"))) +
      guides(colour = guide_legend(override.aes = list(linetype = 0)),
             linetype = guide_legend())
  }  
  plot = plot +
    theme_msstats("VOLCANOPLOT", x.axis.size, y.axis.size,
                  legend.size, strip_background = element_rect(),
                  strip_text_x = element_text(),
                  legend_position = "bottom", legend.title = element_blank())
  plot
}






















#' Preprocess data for comparison plots and create them
#' @inheritParams groupComparisonPlots
#' @param input data.table
#' @param log_base_FC log base for log-fold changes - 2 or 10
#' @keywords internal
.plotComparison = function(
  input, proteins, address, width, height, sig, ylimUp, ylimDown,
  text.angle, dot.size, x.axis.size, y.axis.size, log_base_FC
) {
  adj.pvalue = Protein = ciw = NULL
  
  input = input[!is.na(adj.pvalue), ]
  all_proteins = unique(input$Protein)
  
  if (address == FALSE) {
    if (proteins == "all" | length(proteins) > 1) {
      stop("** Cannnot generate all comparison plots in a screen. Please set one protein at a time.")
    }
  }
  if (proteins != "all") {
    selected_proteins = getSelectedProteins(proteins, all_proteins)
    input = input[Protein %in% selected_proteins, ]
  }
  
  all_proteins = unique(input$Protein)
  input$Protein = factor(input$Protein)
  savePlot(address, "ComparisonPlot", width, height)
  log_fc_column = intersect(colnames(input), c("log2FC", "log10FC"))
  for (i in seq_along(all_proteins)) {
    single_protein = input[Protein == all_proteins[i], ] 		
    single_protein[, ciw := qt(1 - sig / (2 * nrow(single_protein)), single_protein$DF) * single_protein$SE]
    data.table::setnames(single_protein, log_fc_column, "logFC")
    y.limup = ifelse(is.numeric(ylimUp), ylimUp, ceiling(max(single_protein$logFC + single_protein$ciw)))
    y.limdown = ifelse(is.numeric(ylimDown), ylimDown, floor(min(single_protein$logFC - single_protein$ciw)))
    hjust = ifelse(text.angle != 0, 1, 0.5)
    vjust = ifelse(text.angle != 0, 1, 0.5)
    
    plot = .makeComparison(single_protein, log_base_FC, dot.size, x.axis.size,
                           y.axis.size, text.angle, hjust, vjust, y.limdown, 
                           y.limup)
    print(plot)
  }
  if (address != FALSE) {
    dev.off()
  }
}



.makeComparison = function(
  input, log_base, dot.size, x.axis.size, y.axis.size, 
  text.angle, hjust, vjust, y.limdown, y.limup
) {
  logFC = ciw = NULL
  
  protein = unique(input$Protein)
  plot = ggplot(input, aes_string(x = 'Label', y = 'logFC')) +
    geom_errorbar(aes(ymax = logFC + ciw, ymin = logFC - ciw),
                  data = input,
                  width = 0.1,
                  colour = "red") +
    geom_point(size = dot.size, 
               colour = "darkred") +
    scale_x_discrete('Comparison') +
    geom_hline(yintercept = 0, 
               linetype = "twodash", 
               colour = "darkgrey", 
               size = 0.6) +
    labs(title = protein) +
    theme_msstats("COMPARISONPLOT", x.axis.size, y.axis.size, 
                  text_angle = text.angle, text_hjust = hjust, 
                  text_vjust = vjust)
  plot = plot +
    scale_y_continuous(paste0("Log", log_base, "-Fold Change"),
                       limits = c(y.limdown, y.limup))
  plot
}
