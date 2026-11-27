#!/usr/bin/env R
# METABRIC stromal component analysis (deconvolution)
# Matches TCGA deconvolution script structure

packages <- c("dplyr", "ggplot2", "pheatmap", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg); library(pkg, character.only = TRUE)
  }
}

output_dir <- "./metabric_stromal_analysis"
figures_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("METABRIC STROMAL CONTENT ANALYZER\n")
cat("Deconvolution method\n")
cat("========================================\n\n")

# Source all functions from TCGA script with METABRIC-specific modifications
source_tcga_functions <- function() {
  # Reuse TCGA functions: create_reference_profiles, estimate_cell_proportions,
  # categorize_stromal_content, select_high_stromal_tumors, load_target_list,
  # check_gene_overlap, check_targets_in_high_stromal, create_visualizations
  # (Copy those functions here from tcga_deconv.r - they work identically)
}

download_metabric_data <- function() {
  expr_file <- file.path(output_dir, "metabric_expression.rds")
  clinical_file <- file.path(output_dir, "metabric_clinical.rds")
  
  if (file.exists(expr_file)) {
    cat("Loading cached data...\n")
    return(list(expression = readRDS(expr_file), clinical = readRDS(clinical_file)))
  }
  
  cat("Downloading from cBioPortal...\n")
  url <- "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"
  temp_tar <- tempfile(fileext = ".tar.gz")
  temp_dir <- tempdir()
  
  download.file(url, temp_tar, mode = "wb")
  untar(temp_tar, exdir = temp_dir)
  
  # Read expression
  expr_files <- list.files(temp_dir, pattern = "data_mrna.*txt$", full.names = TRUE, recursive = TRUE)
  expr_data <- read.table(expr_files[1], header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  gene_names <- expr_data[[1]]
  expr_data <- expr_data[, -c(1:which(colnames(expr_data) == "Entrez_Gene_Id"))]
  expr_data <- expr_data[!duplicated(gene_names), ]
  rownames(expr_data) <- gene_names[!duplicated(gene_names)]
  
  # Read clinical
  clinical_files <- list.files(temp_dir, pattern = "data_clinical.*txt$", full.names = TRUE, recursive = TRUE)
  clinical_data <- if (length(clinical_files) > 0) {
    read.table(clinical_files[1], header = TRUE, sep = "\t", comment.char = "#", fill = TRUE)
  } else {
    data.frame(SAMPLE_ID = colnames(expr_data))
  }
  
  saveRDS(expr_data, expr_file)
  saveRDS(clinical_data, clinical_file)
  unlink(temp_tar)
  
  return(list(expression = expr_data, clinical = clinical_data))
}

create_reference_profiles <- function() {
  list(
    CAF = list(
      high = c('COL1A1','COL1A2','COL3A1','FAP','ACTA2','PDGFRA','PDGFRB','LUM','DCN','POSTN','THBS2','FN1','VCAN','BGN','SPARC'),
      low = c('PTPRC','CD3D','CD19','EPCAM','KRT8','KRT18','CD34','PECAM1'), weight = 2.0),
    Endothelial = list(
      high = c('PECAM1','CD34','VWF','CDH5','ENG','PLVAP','CLDN5','TEK','FLT1','KDR'),
      low = c('COL1A1','EPCAM','CD3D','FAP','ACTA2'), weight = 1.5),
    Tcells = list(
      high = c('CD3D','CD3E','CD3G','CD8A','CD8B','CD4','IL7R','CD2'),
      low = c('CD19','CD68','EPCAM','COL1A1'), weight = 0.3),
    Bcells = list(
      high = c('CD19','MS4A1','CD79A','CD79B'), low = c('CD3D','CD68','EPCAM'), weight = 0.3),
    Macrophages = list(
      high = c('CD68','CD163','CD14','FCGR3A','ITGAX'), low = c('CD3D','EPCAM'), weight = 0.4),
    NKcells = list(
      high = c('NCAM1','NCR1','KLRB1','NKG7'), low = c('CD3D','CD68'), weight = 0.2),
    Epithelial = list(
      high = c('EPCAM','KRT8','KRT18','KRT19','CDH1','MUC1'),
      low = c('PTPRC','COL1A1','VWF','CD3D'), weight = -1.0)
  )
}

estimate_cell_proportions <- function(expr_data, reference) {
  cat("\n=== Running Pseudo-Deconvolution ===\n")
  cat("Estimating cell type proportions using positive/negative markers\n\n")
  cat("Analyzing", ncol(expr_data), "samples\n\n")
  
  props <- data.frame(sample_id = colnames(expr_data))
  
  for (ct in names(reference)) {
    h <- reference[[ct]]$high[reference[[ct]]$high %in% rownames(expr_data)]
    l <- reference[[ct]]$low[reference[[ct]]$low %in% rownames(expr_data)]
    
    cat(ct, ":\n")
    cat("  High markers found:", length(h), "/", length(reference[[ct]]$high), "\n")
    cat("  Low markers found:", length(l), "/", length(reference[[ct]]$low), "\n")
    
    if (length(h) == 0) { 
      props[[ct]] <- NA
      cat("  WARNING: No markers found!\n\n")
      next 
    }
    
    score <- colMeans(expr_data[h,,drop=F], na.rm=T) - 
             if(length(l)>0) colMeans(expr_data[l,,drop=F], na.rm=T) else 0
    score <- (score - min(score, na.rm=T))/(max(score, na.rm=T) - min(score, na.rm=T))
    props[[ct]] <- score
    
    cat("  Score range:", round(min(score, na.rm=T), 3), "to", round(max(score, na.rm=T), 3), "\n\n")
  }
  
  # Normalize proportions
  cell_cols <- names(reference)[names(reference) %in% names(props)]
  row_sums <- rowSums(props[, cell_cols, drop=F], na.rm=T)
  props[, cell_cols] <- props[, cell_cols] / row_sums
  
  cat("=== Cell Type Abundance (Mean Proportions) ===\n")
  mean_props <- colMeans(props[, cell_cols, drop=F], na.rm=T)
  for (i in seq_along(mean_props)) {
    cat(sprintf("  %-15s: %.3f (%.1f%%)\n", names(mean_props)[i], 
                mean_props[i], mean_props[i] * 100))
  }
  cat("\n")
  
  # Calculate weighted stromal score
  cat("=== Calculating Weighted Stromal Score ===\n")
  props$Stromal_Score <- 0
  for (ct in names(reference)) {
    if (ct %in% names(props) && !all(is.na(props[[ct]]))) {
      props$Stromal_Score <- props$Stromal_Score + props[[ct]] * reference[[ct]]$weight
      cat(sprintf("  %s: weight = %.1f\n", ct, reference[[ct]]$weight))
    }
  }
  
  props$Stromal_Score <- (props$Stromal_Score - min(props$Stromal_Score, na.rm=T))/
                         (max(props$Stromal_Score, na.rm=T) - min(props$Stromal_Score, na.rm=T))
  
  cat("\nStromal Score Statistics:\n")
  cat("  Min:", round(min(props$Stromal_Score, na.rm=T), 3), "\n")
  cat("  Median:", round(median(props$Stromal_Score, na.rm=T), 3), "\n")
  cat("  Mean:", round(mean(props$Stromal_Score, na.rm=T), 3), "\n")
  cat("  Max:", round(max(props$Stromal_Score, na.rm=T), 3), "\n\n")
  
  write.csv(props, file.path(output_dir, "cell_type_proportions.csv"), row.names=F)
  cat("Saved: cell_type_proportions.csv\n")
  
  return(props)
}

categorize_stromal_content <- function(props) {
  cat("\n=== Categorizing Tumors by Stromal Content ===\n")
  
  props <- props[order(props$Stromal_Score, decreasing=T),]
  props$Stromal_Rank <- 1:nrow(props)
  
  q <- quantile(props$Stromal_Score, c(0.33, 0.67), na.rm=T)
  props$Stromal_Category <- cut(props$Stromal_Score, c(-Inf, q, Inf), 
                                labels=c("Low_Stromal","Medium_Stromal","High_Stromal"))
  
  cat("\nStromal Content Categories:\n")
  print(table(props$Stromal_Category))
  
  cat("\nCategory Thresholds:\n")
  cat("  Low:    <", round(q[1], 3), "\n")
  cat("  Medium: ", round(q[1], 3), "-", round(q[2], 3), "\n")
  cat("  High:   >", round(q[2], 3), "\n\n")
  
  write.csv(props, file.path(output_dir, "stromal_content_per_tumor.csv"), row.names=F)
  cat("Saved: stromal_content_per_tumor.csv\n")
  
  return(props)
}

select_high_stromal_tumors <- function(stromal, top_n=NULL, pct=NULL) {
  cat("\n=== Selecting High Stromal Tumors ===\n")
  
  high <- if (!is.null(pct)) {
    threshold <- quantile(stromal$Stromal_Score, pct/100, na.rm=T)
    stromal[stromal$Stromal_Score >= threshold,]
  } else {
    head(stromal, top_n)
  }
  
  if (!is.null(pct)) {
    cat("Selected tumors with stromal score >=", pct, "th percentile\n")
  } else {
    cat("Selected top", nrow(high), "tumors with highest stromal content\n")
  }
  
  cat("\nHigh Stromal Tumors:\n")
  cat("  N samples:", nrow(high), "\n")
  cat("  Mean stromal score:", round(mean(high$Stromal_Score, na.rm=T), 3), "\n")
  cat("  Range:", round(min(high$Stromal_Score, na.rm=T), 3), "-",
      round(max(high$Stromal_Score, na.rm=T), 3), "\n")
  
  # Show cell type composition
  cell_cols <- c("CAF", "Endothelial", "Tcells", "Bcells", "Macrophages", "NKcells", "Epithelial")
  cell_cols <- cell_cols[cell_cols %in% names(high)]
  
  if (length(cell_cols) > 0) {
    cat("\nMean Cell Type Composition in High Stromal Tumors:\n")
    mean_comp <- colMeans(high[, cell_cols, drop=F], na.rm=T)
    for (i in seq_along(mean_comp)) {
      cat(sprintf("  %-15s: %.3f (%.1f%%)\n", names(mean_comp)[i], 
                  mean_comp[i], mean_comp[i] * 100))
    }
  }
  
  write.csv(high, file.path(output_dir, "high_stromal_tumors.csv"), row.names=F)
  cat("\nSaved: high_stromal_tumors.csv\n")
  
  return(high)
}

load_target_list <- function(file=NULL) {
  cat("\n=== Loading Target Gene List (TEST GENES) ===\n")
  
  if (is.null(file) || file == "") {
    cat("No file provided. Using example targets...\n")
    targets <- c('COL1A1','FAP','EPCAM','KRT8','CD3D')
    writeLines(targets, file.path(output_dir, "example_targets.txt"))
    cat("Loaded", length(targets), "example genes\n")
    return(targets)
  }
  
  if (!file.exists(file)) {
    cat("\nERROR: File not found! Using example targets...\n")
    return(load_target_list(NULL))
  }
  
  cat("Loading from:", file, "\n")
  
  targets <- if (grepl("\\.(csv|tsv)$", file)) {
    as.character(read.csv(file, stringsAsFactors=FALSE)[,1])
  } else {
    # Use readLines with warn=FALSE to suppress incomplete line warning
    readLines(file, warn=FALSE)
  }
  
  targets <- trimws(targets)
  targets <- targets[targets != ""]
  targets <- targets[!grepl("^#", targets)]
  targets <- unique(targets)
  
  cat("Loaded", length(targets), "target genes\n")
  
  return(targets)
}

check_gene_overlap <- function(ref, tgts) {
  cat("\n=== Checking for Gene Overlap ===\n")
  
  ref_genes <- unique(unlist(lapply(ref, function(x) c(x$high, x$low))))
  overlapping <- tgts[tgts %in% ref_genes]
  
  if (length(overlapping) > 0) {
    cat("\n*** WARNING: OVERLAP DETECTED ***\n")
    cat("These test genes are in the reference signatures:\n")
    for (gene in overlapping) {
      cat("  -", gene, "\n")
    }
    cat("\nThis may cause circular analysis!\n\n")
  } else {
    cat("\nNo overlap detected - analysis is independent\n\n")
  }
  
  return(overlapping)
}

check_targets_in_high_stromal <- function(expr, tgts, high, stromal) {
  cat("\n=== Testing Target Genes in High vs Low Stromal Tumors ===\n")
  
  avail <- tgts[tgts %in% rownames(expr)]
  missing <- tgts[!tgts %in% rownames(expr)]
  
  cat("\nTarget Gene Status:\n")
  cat("  Found:", length(avail), "/", length(tgts), "\n")
  cat("  Missing:", length(missing), "\n")
  
  if (length(missing) > 0 && length(missing) <= 20) {
    cat("\nMissing targets:\n")
    for (gene in missing) cat("  ", gene, "\n")
  }
  
  if (length(avail) == 0) {
    cat("\nNo target genes found in dataset!\n")
    return(NULL)
  }
  
  h_samp <- high$sample_id
  l_samp <- tail(stromal, nrow(high))$sample_id
  
  cat("\n=== Comparison Groups ===\n")
  cat("High stromal tumors:", length(h_samp), "\n")
  cat("Low stromal tumors:", length(l_samp), "\n\n")
  
  stats <- data.frame(
    gene = avail,
    high_stromal_mean = rowMeans(expr[avail, h_samp, drop=F], na.rm=T),
    high_stromal_median = apply(expr[avail, h_samp, drop=F], 1, median, na.rm=T),
    low_stromal_mean = rowMeans(expr[avail, l_samp, drop=F], na.rm=T),
    low_stromal_median = apply(expr[avail, l_samp, drop=F], 1, median, na.rm=T),
    fold_change = NA,
    p_value = NA,
    effect_size = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(stats)) {
    h <- expr[stats$gene[i], h_samp]; l <- expr[stats$gene[i], l_samp]
    h <- h[!is.na(h)]; l <- l[!is.na(l)]
    
    stats$fold_change[i] <- log2((stats$high_stromal_mean[i] + 0.01)/(stats$low_stromal_mean[i] + 0.01))
    
    pooled_sd <- sqrt((sd(h)^2 + sd(l)^2)/2)
    stats$effect_size[i] <- if (pooled_sd > 0) (mean(h) - mean(l))/pooled_sd else 0
    
    test <- tryCatch(wilcox.test(h, l), error=function(e) NULL)
    if (!is.null(test)) stats$p_value[i] <- test$p.value
  }
  
  stats$p_adjusted <- p.adjust(stats$p_value, "BH")
  stats$significant <- ifelse(stats$p_adjusted < 0.05 & abs(stats$effect_size) > 0.5, "***",
                             ifelse(stats$p_adjusted < 0.05 & abs(stats$effect_size) > 0.2, "**",
                                   ifelse(stats$p_adjusted < 0.05, "*", "ns")))
  
  stats <- stats[order(stats$fold_change, decreasing=T),]
  
  cat("Target Expression in High vs Low Stromal Tumors:\n")
  cat(strrep("-", 110), "\n")
  cat(sprintf("%-12s %10s %10s %8s %10s %10s %8s %5s\n",
              "Gene", "High Mean", "Low Mean", "Log2FC", "P-adj", "Effect", "Sig", ""))
  cat(strrep("-", 110), "\n")
  
  for (i in 1:nrow(stats)) {
    cat(sprintf("%-12s %10.2f %10.2f %8.2f %10.2e %10.2f %8s %5s\n",
                stats$gene[i],
                stats$high_stromal_mean[i],
                stats$low_stromal_mean[i],
                stats$fold_change[i],
                stats$p_adjusted[i],
                stats$effect_size[i],
                ifelse(abs(stats$effect_size[i]) > 0.5, "LARGE",
                       ifelse(abs(stats$effect_size[i]) > 0.2, "SMALL", "TINY")),
                stats$significant[i]))
  }
  
  cat("\nInterpretation Guide:\n")
  cat("  Effect Size: >0.8 = Large, 0.5-0.8 = Medium, 0.2-0.5 = Small, <0.2 = Negligible\n")
  cat("  *** = p<0.05 + large effect, ** = p<0.05 + small effect, * = p<0.05 only\n\n")
  
  write.csv(stats, file.path(output_dir, "targets_in_high_stromal_tumors.csv"), row.names=F)
  cat("Saved: targets_in_high_stromal_tumors.csv\n")
  
  return(list(stats=stats, high_stromal_expr=expr[avail, h_samp, drop=F], 
             low_stromal_expr=expr[avail, l_samp, drop=F]))
}

create_visualizations <- function(stromal, high, tgt_result) {
  png(file.path(figures_dir, "distribution.png"), w=12, h=6, units="in", res=300)
  par(mfrow=c(1,2))
  hist(stromal$Stromal_Score, breaks=30, col="#3498db", main="Stromal Score", xlab="Score")
  barplot(table(stromal$Stromal_Category), col=c("#2ecc71","#f39c12","#e74c3c"), main="Categories")
  dev.off()
  
  if (!is.null(tgt_result)) {
    png(file.path(figures_dir, "targets.png"), w=10, h=max(6, nrow(tgt_result$stats)*0.3), units="in", res=300)
    par(mar=c(5,10,4,2))
    s <- tgt_result$stats[order(tgt_result$stats$fold_change),]
    cols <- ifelse(s$p_adjusted < 0.05 & s$fold_change > 0 & abs(s$effect_size) > 0.5, "#e74c3c",
                  ifelse(s$p_adjusted < 0.05 & s$fold_change < 0 & abs(s$effect_size) > 0.5, "#2ecc71", "#95a5a6"))
    barplot(s$fold_change, names=s$gene, horiz=T, las=1, col=cols, main="Target Expression: High vs Low")
    abline(v=0, lty=2)
    dev.off()
  }
}

# Command line: Rscript metabric_deconv.R <gene list> [<high-stromal selection>]
# Without arguments the script asks for these instead.
cli_args <- commandArgs(trailingOnly = TRUE)
answer <- function(i) {
  if (length(cli_args) == 0) return(readline())
  value <- if (length(cli_args) >= i) cli_args[i] else ""
  cat(value, "\n")
  value
}

main <- function() {
  metabric <- download_metabric_data()
  expr <- metabric$expression
  reference <- create_reference_profiles()
  props <- estimate_cell_proportions(expr, reference)
  stromal <- categorize_stromal_content(props)
  
  cat("\n", strrep("=", 70), "\n")
  cat("How many high stromal tumors to select?\n")
  cat("Enter number (e.g., 100) or percentile (e.g., 75 for top 25%): ")
  sel <- suppressWarnings(as.numeric(answer(2)))
  if (is.na(sel)) sel <- 100
  
  high <- if (sel > 50) select_high_stromal_tumors(stromal, pct=sel) 
          else select_high_stromal_tumors(stromal, top_n=sel)
  
  cat("\n", strrep("=", 70), "\n")
  cat("Enter path to your target gene list file\n")
  cat("(or press Enter to use example targets): ")
  tfile <- answer(1)
  tfile <- trimws(gsub('^["\']|["\']$', '', tfile))
  
  targets <- load_target_list(tfile)
  
  overlapping <- check_gene_overlap(reference, targets)
  if (length(overlapping) > 0) {
    cat("\n", strrep("=", 70), "\n")
    cat("WARNING: Some test genes overlap with reference signatures!\n")
    cat("Continue despite overlap? (y/n): ")
    if (tolower(trimws(readline())) != "y") stop("Analysis canceled due to gene overlap")
  }
  
  tgt_result <- check_targets_in_high_stromal(expr, targets, high, stromal)
  create_visualizations(stromal, high, tgt_result)
  
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE!\n")
  cat(strrep("=", 70), "\n")
  cat("\nOutput directory:", normalizePath(output_dir), "\n")
  cat("\nFiles created:\n")
  cat("  cell_type_proportions.csv\n")
  cat("  stromal_content_per_tumor.csv\n")
  cat("  high_stromal_tumors.csv\n")
  cat("  targets_in_high_stromal_tumors.csv\n")
  cat("  figures/distribution.png\n")
  cat("  figures/targets.png\n\n")
  
  cat("DECONVOLUTION ANALYSIS NOTES:\n")
  cat("  - Used pseudo-deconvolution with positive/negative markers\n")
  cat("  - Cell type proportions estimated for: CAF, Endothelial, Immune cells\n")
  cat("  - Stromal score = weighted sum (CAF×2.0 + Endothelial×1.5 - Epithelial×1.0)\n")
  cat("  - Higher weight on CAF reflects its dominant stromal role\n")
  cat("  - Test genes analyzed independently from reference signatures\n\n")
  
  cat("INTERPRETATION:\n")
  cat("  - High stromal tumors = high CAF + endothelial, low epithelial\n")
  cat("  - Low stromal tumors = high epithelial, low CAF + endothelial\n")
  cat("  - Effect sizes show biological significance (>0.5 = meaningful)\n\n")
}

main()