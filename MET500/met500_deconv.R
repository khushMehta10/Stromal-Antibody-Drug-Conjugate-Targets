#!/usr/bin/env R
# MET500 BRCA stromal component analysis (deconvolution)
# Matches TCGA deconvolution script with interactive interface

packages <- c("dplyr", "ggplot2", "pheatmap", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

output_dir <- "./MET500_stromal_analysis"
figures_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("MET500 BRCA STROMAL CONTENT ANALYZER\n")
cat("Deconvolution method\n")
cat("========================================\n\n")

################################################################################
# Load Data with Interactive Prompt
################################################################################

load_met500_data <- function() {
  cat("\n=== Loading MET500 Expression Data ===\n")
  cat("\nEnter path to M.mx.log2.zscore.txt file: ")
  expr_file <- trimws(gsub('^["\']|["\']$', '', answer(3)))
  
  if (expr_file == "") expr_file <- "M.mx.log2.zscore.txt"
  if (!file.exists(expr_file)) stop("File not found: ", expr_file)
  
  cat("Loading...\n")
  expr_data <- read.table(expr_file, header = TRUE, row.names = 1, 
                         sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Clean ENSEMBL IDs
  if (any(grepl("ENSG", rownames(expr_data)))) {
    rownames(expr_data) <- sub("\\..*", "", rownames(expr_data))
    expr_data <- expr_data[!duplicated(rownames(expr_data)), ]
  }
  
  cat("Loaded:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
  saveRDS(as.matrix(expr_data), file.path(output_dir, "expression.rds"))
  return(as.matrix(expr_data))
}

################################################################################
# Reference Profiles (ENSEMBL IDs)
################################################################################

create_reference_profiles <- function() {
  list(
    CAF = list(
      high = c('ENSG00000108821','ENSG00000164692','ENSG00000078098','ENSG00000134853','ENSG00000196712'),
      low = c('ENSG00000081237','ENSG00000167286','ENSG00000119888'), weight = 2.0),
    Endothelial = list(
      high = c('ENSG00000261371','ENSG00000174059','ENSG00000110799','ENSG00000179776'),
      low = c('ENSG00000108821','ENSG00000119888'), weight = 1.5),
    Tcells = list(
      high = c('ENSG00000167286','ENSG00000198851','ENSG00000153563','ENSG00000010610'),
      low = c('ENSG00000177455','ENSG00000091972'), weight = 0.3),
    Epithelial = list(
      high = c('ENSG00000119888','ENSG00000170421','ENSG00000111057','ENSG00000039068'),
      low = c('ENSG00000081237','ENSG00000108821'), weight = -1.0)
  )
}

################################################################################
# Deconvolution
################################################################################

estimate_cell_proportions <- function(expr_data, reference) {
  cat("\n=== Running Pseudo-Deconvolution ===\n")
  props <- data.frame(sample_id = colnames(expr_data))
  
  for (ct in names(reference)) {
    h <- reference[[ct]]$high[reference[[ct]]$high %in% rownames(expr_data)]
    l <- reference[[ct]]$low[reference[[ct]]$low %in% rownames(expr_data)]
    
    if (length(h) == 0) { props[[ct]] <- NA; next }
    
    score <- colMeans(expr_data[h,,drop=F]) - 
             if(length(l)>0) colMeans(expr_data[l,,drop=F]) else 0
    props[[ct]] <- (score - min(score))/(max(score) - min(score))
  }
  
  # Weighted stromal score
  props$Stromal_Score <- 0
  for (ct in names(reference)) {
    if (ct %in% names(props) && !all(is.na(props[[ct]]))) {
      props$Stromal_Score <- props$Stromal_Score + props[[ct]] * reference[[ct]]$weight
    }
  }
  props$Stromal_Score <- (props$Stromal_Score - min(props$Stromal_Score))/
                         (max(props$Stromal_Score) - min(props$Stromal_Score))
  
  write.csv(props, file.path(output_dir, "cell_proportions.csv"), row.names=F)
  return(props)
}

################################################################################
# Categorize & Select
################################################################################

categorize_stromal <- function(props) {
  props <- props[order(props$Stromal_Score, decreasing=T),]
  q <- quantile(props$Stromal_Score, c(0.33, 0.67))
  props$Category <- cut(props$Stromal_Score, c(-Inf, q, Inf), 
                       labels=c("Low","Medium","High"))
  write.csv(props, file.path(output_dir, "stromal_categories.csv"), row.names=F)
  return(props)
}

select_high_stromal <- function(stromal, top_n=NULL, pct=NULL) {
  if (!is.null(pct)) {
    thresh <- quantile(stromal$Stromal_Score, pct/100)
    high <- stromal[stromal$Stromal_Score >= thresh,]
  } else {
    high <- head(stromal, top_n)
  }
  write.csv(high, file.path(output_dir, "high_stromal.csv"), row.names=F)
  return(high)
}

################################################################################
# Target Analysis
################################################################################

analyze_targets <- function(expr_data, targets, high, stromal, mapping=NULL) {
  avail <- targets[targets %in% rownames(expr_data)]
  if (length(avail)==0) return(NULL)
  
  high_samples <- high$sample_id
  low_samples <- tail(stromal, nrow(high))$sample_id
  
  stats <- data.frame(
    ENSEMBL_ID = avail,
    Gene = avail,
    High_Mean = rowMeans(expr_data[avail, high_samples, drop=F]),
    Low_Mean = rowMeans(expr_data[avail, low_samples, drop=F]),
    Diff = NA, P_val = NA, Effect = NA
  )
  
  # Map symbols
  if (!is.null(mapping)) {
    for (i in 1:nrow(stats)) {
      m <- match(stats$ENSEMBL_ID[i], mapping$ENSEMBL_ID)
      if (!is.na(m)) stats$Gene[i] <- mapping$Gene_Symbol[m]
    }
  }
  
  for (i in 1:nrow(stats)) {
    h <- expr_data[stats$ENSEMBL_ID[i], high_samples]
    l <- expr_data[stats$ENSEMBL_ID[i], low_samples]
    stats$Diff[i] <- mean(h) - mean(l)
    stats$Effect[i] <- (mean(h)-mean(l))/sqrt((sd(h)^2+sd(l)^2)/2)
    test <- tryCatch(wilcox.test(h,l), error=function(e) NULL)
    if (!is.null(test)) stats$P_val[i] <- test$p.value
  }
  
  stats$P_adj <- p.adjust(stats$P_val, "BH")
  stats$Sig <- ifelse(stats$P_adj<0.05 & abs(stats$Effect)>0.5, "***",
                     ifelse(stats$P_adj<0.05 & abs(stats$Effect)>0.2, "**",
                           ifelse(stats$P_adj<0.05, "*", "ns")))
  
  stats <- stats[order(stats$Effect, decreasing=T),]
  write.csv(stats, file.path(output_dir, "target_analysis.csv"), row.names=F)
  
  cat("\nTarget Analysis Results:\n")
  print(stats[,c("Gene","Diff","Effect","P_adj","Sig")])
  
  return(list(stats=stats, 
             high_expr=expr_data[avail, high_samples, drop=F],
             low_expr=expr_data[avail, low_samples, drop=F]))
}

################################################################################
# Visualizations
################################################################################

create_visualizations <- function(stromal, high, target_result) {
  # Distribution plot
  png(file.path(figures_dir, "stromal_distribution.png"), w=14, h=6, units="in", res=300)
  par(mfrow=c(1,2))
  hist(stromal$Stromal_Score, breaks=30, col="#3498db", main="Stromal Score Distribution")
  barplot(table(stromal$Category), col=c("#2ecc71","#f39c12","#e74c3c"), 
         main="Categories")
  dev.off()
  
  # Target fold change
  if (!is.null(target_result)) {
    stats <- target_result$stats
    png(file.path(figures_dir, "target_foldchange.png"), w=10, h=max(6,nrow(stats)*0.3), 
       units="in", res=300)
    par(mar=c(5,10,4,2))
    cols <- ifelse(stats$P_adj<0.05 & stats$Diff>0 & abs(stats$Effect)>0.5, "#e74c3c",
                  ifelse(stats$P_adj<0.05 & stats$Diff<0 & abs(stats$Effect)>0.5, "#2ecc71", "#95a5a6"))
    barplot(stats$Diff, names=stats$Gene, horiz=T, las=1, col=cols, 
           main="Target Expression: High vs Low Stromal")
    abline(v=0, lty=2)
    dev.off()
  }
}

################################################################################
# Main
################################################################################

# Command line: Rscript met500_deconv.R <gene list> [<high-stromal selection>] [<expression file>]
# Without arguments the script asks for these instead.
cli_args <- commandArgs(trailingOnly = TRUE)
answer <- function(i) {
  if (length(cli_args) == 0) return(readline())
  value <- if (length(cli_args) >= i) cli_args[i] else ""
  cat(value, "\n")
  value
}

main <- function() {
  expr_data <- load_met500_data()
  reference <- create_reference_profiles()
  props <- estimate_cell_proportions(expr_data, reference)
  stromal <- categorize_stromal(props)
  
  cat("\nSelect high stromal samples (enter number or percentile >50): ")
  sel <- suppressWarnings(as.numeric(answer(2)))
  if (is.na(sel)) sel <- 100
  high <- if (sel > 50) select_high_stromal(stromal, pct=sel) 
          else select_high_stromal(stromal, top_n=sel)
  
  cat("\nEnter target gene file path (or press Enter for defaults): ")
  tfile <- trimws(gsub('^["\']|["\']$', '', answer(1)))
  
  targets <- if (tfile == "") {
    c('ENSG00000108821','ENSG00000078098','ENSG00000119888','ENSG00000170421')
  } else {
    readLines(tfile)
  }
  
  mapping <- data.frame(
    ENSEMBL_ID=c('ENSG00000108821','ENSG00000078098','ENSG00000119888','ENSG00000170421'),
    Gene_Symbol=c('COL1A1','FAP','EPCAM','KRT8')
  )
  
  target_result <- analyze_targets(expr_data, targets, high, stromal, mapping)
  create_visualizations(stromal, high, target_result)
  
  cat("\n========================================\n")
  cat("ANALYSIS COMPLETE!\n")
  cat("Files saved in:", normalizePath(output_dir), "\n")
  cat("========================================\n\n")
}

main()