## Original script by Thomas Kuilman
## Modified for fgsea by Christina B. Pedersen

replotFGSEA <- function(stats, gene.set.name, gene.sets, fgsea, class.name, metric.range = NULL,
                       enrichment.score.range = NULL) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
  statsAdj <- statsAdj/max(abs(statsAdj))
  
  pathway <- gene.sets[[gene.set.name]]
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(-30, xs, n + 31), y = c(0, ys, 0))
  
  if (is.null(enrichment.score.range)) {
    enrichment.score.range <- c(min(toPlot$y), max(toPlot$y)) 
  }
  metric.range <- c(min(stats[ord]), max(stats[ord]))
  
  # Get fgsea results for gene set
  fgsea.set <- fgsea[which(fgsea$pathway==gene.set.name),]
  
  # Get enrichment score
  gsea.enrichment.score <- round(fgsea.set$ES,3)

  # Get gene set name
  gsea.normalized.enrichment.score <- round(fgsea.set$NES,3)

  # Get nominal p-value
  gsea.p.value <- round(fgsea.set$pval,3)

  # Get FDR
  gsea.fdr <- round(fgsea.set$padj,3)

  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  #dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  #layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  
  plot(toPlot$x,
       toPlot$y, type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list(gene.set.name, font = 1, cex = 1.5),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(n * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR adj. p-value:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0))
  } else {
    text(n * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR adj. p-value:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(-30, n+31))
  abline(v = pathway, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- stats[ord] - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(-30, n+31))
  box()
  text(n / 2, 0.7,
       labels = class.name)
  text(n * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(n * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(stats[ord])
  plot(stats[ord], type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(-30, n+31),
       ylim = metric.range, yaxs = "i",
       ylab = "Ranking metric", 
       panel.first = abline(h = seq(metric.range[1] / 2,
                                    metric.range[2] - metric.range[1] / 4,
                                    metric.range[2] / 2), col = "gray95", lty = 2))
  
  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlim = c(-30, n+31),  xaxt='n',
          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
          ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
  
  abline(v = length(which(rank.metric$values > 0)), lty = 2, col = 'grey20')
  text(x = length(which(rank.metric$values > 0))*1.25, y = max(rank.metric$values)*0.2, labels = paste0('Zero cross at ', format(length(which(rank.metric$values > 0)), nsmall=0, big.mark=",")))
  
  box()
  
  # Reset to default
  par(def.par)
  
}



replot_multiFGSEA <- function(stats, gene.sets, fgsea, class.name, metric.range = NULL,
                       enrichment.score.range = NULL) {
  
  # Processing genes (global for dataset)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
  statsAdj <- statsAdj/max(abs(statsAdj))
  
  metric.range <- c(min(stats[ord]), max(stats[ord]))
  
  # Processing each gene set and saving results
  toPlot <- pathways <- gsea.enrichment.score <- list()
  for (gene.set.name in names(gene.sets)) {
    pathway <- gene.sets[[gene.set.name]]
    
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathways[[gene.set.name]] <- pathway <- sort(pathway)
    
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot[[gene.set.name]] <- data.frame(x = c(-30, xs, n + 31), y = c(0, ys, 0))
    
    # Get fgsea results for gene set
    fgsea.set <- fgsea[which(fgsea$pathway==gene.set.name),]
    
    # Get enrichment score
    gsea.enrichment.score[[gene.set.name]] <- round(fgsea.set$ES,3)
    
    # # Get gene set name
    # gsea.normalized.enrichment.score <- round(fgsea.set$NES,3)
    # 
    # # Get nominal p-value
    # gsea.p.value <- round(fgsea.set$pval,3)
    # 
    # # Get FDR
    # gsea.fdr <- round(fgsea.set$padj,3)
  }
  
  # Range for scores
  if (is.null(enrichment.score.range)) {
    enrichment.score.range <- c(min(sapply(toPlot, function(d) {min(d[,'y'])})), max(sapply(toPlot, function(d) {max(d[,'y'])}))) 
  }
  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  #dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(seq(1:(3+length(gene.sets)))), heights = c(1.7, rep(0.3, length(gene.sets)), 0.2, 2))
  #layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  
  colors <- RColorBrewer::brewer.pal(length(toPlot), 'Set1')
  plot(toPlot[[1]]$x,
       toPlot[[1]]$y, type = "l", col = colors[1], lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list('GSEA enrichment plot', font = 1, cex = 1.5),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  for (i in 2:length(toPlot)) {
    lines(toPlot[[i]]$x,
          toPlot[[i]]$y, type = 'l', col = colors[i], lwd = 1.5)
  }
  
  plot.coordinates <- par("usr")
  text(n * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
       paste("Enrinchment Scores:\n",
             paste(names(gsea.enrichment.score), gsea.enrichment.score, collapse = '\n')), cex = 0.9, adj = c(1, 1))
  
  # Barcodes
  for (i in 1:length(pathways)) {
    par(mar = c(0, 5, 0, 2))
    plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(-30, n+31))
    abline(v = pathways[[i]], lwd = 0.75, col = colors[i])
  }
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- stats[ord] - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(-30, n+31))
  box()
  text(n / 2, 0.7,
       labels = class.name)
  text(n * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(n * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(stats[ord])
  plot(stats[ord], type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(-30, n+31),
       ylim = metric.range, yaxs = "i",
       ylab = "Ranking metric", 
       panel.first = abline(h = seq(metric.range[1] / 2,
                                    metric.range[2] - metric.range[1] / 4,
                                    metric.range[2] / 2), col = "gray95", lty = 2))
  
  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlim = c(-30, n+31),  xaxt='n',
          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
          ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
  
  abline(v = length(which(rank.metric$values > 0)), lty = 2, col = 'grey20')
  text(x = length(which(rank.metric$values > 0))*1.3, y = max(rank.metric$values)*0.2, labels = paste0('Zero cross at ', format(length(which(rank.metric$values > 0)), nsmall=0, big.mark=",")))
  
  box()
  
  # Reset to default
  par(def.par)
  
}


