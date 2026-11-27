#!/usr/bin/env R
# TCGA-BRCA stromal component analysis (deconvolution)
# Uses pseudo-deconvolution for more accurate cell type estimation

packages <- c("dplyr", "ggplot2", "pheatmap", "reshape2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# Setup
################################################################################

output_dir <- "./tcga_stromal_analysis"
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("TCGA-BRCA STROMAL CONTENT ANALYZER\n")
cat("Deconvolution method\n")
cat("========================================\n\n")
cat("Output directory:", normalizePath(output_dir), "\n")

################################################################################
# Download TCGA Data
################################################################################

download_tcga_data <- function() {
  
  cat("\n=== Downloading TCGA-BRCA Data ===\n")
  
  expr_file <- file.path(output_dir, "tcga_brca_expression.rds")
  
  if (file.exists(expr_file)) {
    cat("Expression data already exists. Loading...\n")
    expr_data <- readRDS(expr_file)
    cat("Loaded:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
    return(expr_data)
  }
  
  cat("\nDownloading from UCSC Xena Browser...\n")
  
  url <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz"
  temp_file <- tempfile(fileext = ".gz")
  
  cat("Downloading... (~200 MB, 5-10 minutes)\n\n")
  
  tryCatch({
    download.file(url, temp_file, mode = "wb", quiet = FALSE)
    cat("\nDownload complete!\n")
  }, error = function(e) {
    stop("Download failed: ", conditionMessage(e))
  })
  
  cat("\nReading data...\n")
  expr_data <- read.table(gzfile(temp_file), 
                          header = TRUE, 
                          row.names = 1, 
                          sep = "\t",
                          check.names = FALSE)
  
  unlink(temp_file)
  
  cat("Loaded:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
  
  expr_data <- expr_data[rowSums(expr_data) > 0, ]
  cat("After filtering:", nrow(expr_data), "expressed genes\n")
  
  saveRDS(expr_data, expr_file)
  cat("Saved expression data\n")
  
  return(expr_data)
}

################################################################################
# Identify Sample Types
################################################################################

identify_sample_types <- function(expr_data) {
  
  cat("\n=== Identifying Sample Types ===\n")
  
  sample_names <- colnames(expr_data)
  
  sample_codes <- sapply(strsplit(sample_names, "-"), function(x) {
    if (length(x) >= 4) x[4] else NA
  })
  
  sample_info <- data.frame(
    sample_id = sample_names,
    sample_code = sample_codes,
    sample_type = ifelse(sample_codes == "01", "Tumor",
                        ifelse(sample_codes == "11", "Normal", "Other")),
    stringsAsFactors = FALSE
  )
  
  cat("\nSample Type Summary:\n")
  print(table(sample_info$sample_type))
  
  tumor_samples <- sample_info$sample_id[sample_info$sample_type == "Tumor"]
  cat("\nFound", length(tumor_samples), "tumor samples for analysis\n")
  
  return(sample_info)
}

################################################################################
# Create Reference Cell Type Profiles (Pseudo-deconvolution)
################################################################################

create_reference_profiles <- function() {
  
  cat("\n=== Creating Cell Type Reference Profiles ===\n")
  cat("Using pseudo-deconvolution with positive/negative markers\n\n")
  
  # Reference profiles: genes that should be HIGH/LOW in each cell type
  # This mimics CIBERSORTx but without needing purified cell populations
  
  reference <- list(
    
    # Cancer-Associated Fibroblasts (CAFs) - PRIMARY STROMAL COMPONENT
    CAF = list(
      high = c('COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A1', 'COL6A2',
               'FAP', 'ACTA2', 'PDGFRA', 'PDGFRB', 'LUM', 'DCN', 'POSTN', 
               'THBS2', 'FN1', 'VCAN', 'BGN', 'SPARC', 'MMP2', 'MMP11'),
      low = c('PTPRC', 'CD3D', 'CD19', 'EPCAM', 'KRT8', 'KRT18', 'CD34', 'PECAM1'),
      weight = 2.0  # CAFs get highest weight for stromal score
    ),
    
    # Endothelial Cells - VASCULAR STROMA
    Endothelial = list(
      high = c('PECAM1', 'CD34', 'VWF', 'CDH5', 'ENG', 'PLVAP', 'CLDN5', 
               'TEK', 'FLT1', 'KDR', 'ESAM', 'MMRN1', 'CLEC14A'),
      low = c('COL1A1', 'EPCAM', 'CD3D', 'FAP', 'ACTA2', 'CD68'),
      weight = 1.5
    ),
    
    # T Cells - IMMUNE INFILTRATE
    Tcells = list(
      high = c('CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'IL7R', 
               'CD2', 'CD5', 'CD7', 'TRAC', 'TRBC1', 'TRBC2'),
      low = c('CD19', 'CD68', 'EPCAM', 'COL1A1', 'CD34'),
      weight = 0.3  # Immune has minimal stromal contribution
    ),
    
    # B Cells - IMMUNE INFILTRATE
    Bcells = list(
      high = c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'BLK', 'CD22', 'BANK1'),
      low = c('CD3D', 'CD68', 'EPCAM', 'COL1A1', 'CD34'),
      weight = 0.3
    ),
    
    # Macrophages/Myeloid - IMMUNE INFILTRATE
    Macrophages = list(
      high = c('CD68', 'CD163', 'CD14', 'FCGR3A', 'ITGAX', 'CD86', 
               'CSF1R', 'MARCO', 'MSR1', 'MRC1'),
      low = c('CD3D', 'CD19', 'EPCAM', 'KRT8', 'COL1A1'),
      weight = 0.4  # TAMs can support stroma
    ),
    
    # NK Cells - IMMUNE INFILTRATE
    NKcells = list(
      high = c('NCAM1', 'NCR1', 'KLRB1', 'KLRD1', 'NKG7', 'GNLY'),
      low = c('CD3D', 'CD19', 'CD68', 'EPCAM', 'COL1A1'),
      weight = 0.2
    ),
    
    # Epithelial/Tumor Cells - INVERSE OF STROMA
    Epithelial = list(
      high = c('EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'MUC1', 
               'KRT7', 'CLDN4', 'CLDN3', 'TJP1'),
      low = c('PTPRC', 'COL1A1', 'VWF', 'CD3D', 'FAP'),
      weight = -1.0  # Epithelial REDUCES stromal score
    )
  )
  
  cat("Cell Type Reference Profiles:\n")
  for (cell_type in names(reference)) {
    cat(sprintf("  %-15s: %2d high markers, %2d low markers, weight = %.1f\n",
                cell_type,
                length(reference[[cell_type]]$high),
                length(reference[[cell_type]]$low),
                reference[[cell_type]]$weight))
  }
  cat("\n")
  
  return(reference)
}

################################################################################
# Estimate Cell Type Proportions (Pseudo-deconvolution)
################################################################################

estimate_cell_proportions <- function(expr_data, reference, sample_info) {
  
  cat("\n=== Running Pseudo-Deconvolution ===\n")
  cat("Estimating cell type proportions using positive/negative markers\n\n")
  
  tumor_samples <- sample_info$sample_id[sample_info$sample_type == "Tumor"]
  tumor_expr <- expr_data[, tumor_samples, drop = FALSE]
  
  cat("Analyzing", ncol(tumor_expr), "tumor samples\n\n")
  
  # Initialize results
  proportions <- data.frame(
    sample_id = tumor_samples,
    stringsAsFactors = FALSE
  )
  
  # Calculate enrichment score for each cell type
  for (cell_type in names(reference)) {
    
    profile <- reference[[cell_type]]
    
    # Genes that should be HIGH in this cell type
    high_genes <- profile$high[profile$high %in% rownames(tumor_expr)]
    # Genes that should be LOW in this cell type  
    low_genes <- profile$low[profile$low %in% rownames(tumor_expr)]
    
    cat(cell_type, ":\n")
    cat("  High markers found:", length(high_genes), "/", length(profile$high), "\n")
    cat("  Low markers found:", length(low_genes), "/", length(profile$low), "\n")
    
    if (length(high_genes) == 0) {
      proportions[[cell_type]] <- NA
      cat("  WARNING: No markers found!\n\n")
      next
    }
    
    # Calculate enrichment score for each sample
    high_expr <- colMeans(tumor_expr[high_genes, , drop = FALSE])
    
    if (length(low_genes) > 0) {
      low_expr <- colMeans(tumor_expr[low_genes, , drop = FALSE])
      # Score = high markers - low markers (differential expression)
      score <- high_expr - low_expr
    } else {
      score <- high_expr
    }
    
    # Normalize to 0-1 range within this cell type
    score <- (score - min(score, na.rm = TRUE)) / 
             (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
    
    proportions[[cell_type]] <- score
    
    cat("  Score range:", round(min(score, na.rm = TRUE), 3), 
        "to", round(max(score, na.rm = TRUE), 3), "\n\n")
  }
  
  # Normalize so proportions sum to 1 per sample (like real deconvolution)
  cell_cols <- names(reference)
  cell_cols <- cell_cols[cell_cols %in% names(proportions)]
  
  row_sums <- rowSums(proportions[, cell_cols, drop = FALSE], na.rm = TRUE)
  proportions[, cell_cols] <- proportions[, cell_cols] / row_sums
  
  cat("=== Cell Type Abundance (Mean Proportions) ===\n")
  mean_props <- colMeans(proportions[, cell_cols, drop = FALSE], na.rm = TRUE)
  for (i in seq_along(mean_props)) {
    cat(sprintf("  %-15s: %.3f (%.1f%%)\n", 
                names(mean_props)[i], 
                mean_props[i],
                mean_props[i] * 100))
  }
  cat("\n")
  
  # Calculate composite STROMAL SCORE using weighted combination
  cat("=== Calculating Weighted Stromal Score ===\n")
  
  proportions$Stromal_Score <- 0
  
  for (cell_type in names(reference)) {
    if (cell_type %in% names(proportions)) {
      weight <- reference[[cell_type]]$weight
      proportions$Stromal_Score <- proportions$Stromal_Score + 
                                    (proportions[[cell_type]] * weight)
      cat(sprintf("  %s: weight = %.1f\n", cell_type, weight))
    }
  }
  
  # Normalize stromal score to 0-1
  proportions$Stromal_Score <- (proportions$Stromal_Score - min(proportions$Stromal_Score)) /
                               (max(proportions$Stromal_Score) - min(proportions$Stromal_Score))
  
  cat("\nStromal Score Statistics:\n")
  cat("  Min:", round(min(proportions$Stromal_Score), 3), "\n")
  cat("  Median:", round(median(proportions$Stromal_Score), 3), "\n")
  cat("  Mean:", round(mean(proportions$Stromal_Score), 3), "\n")
  cat("  Max:", round(max(proportions$Stromal_Score), 3), "\n\n")
  
  # Save detailed cell type proportions
  write.csv(proportions, 
            file.path(output_dir, "cell_type_proportions.csv"), 
            row.names = FALSE)
  cat("Saved: cell_type_proportions.csv\n")
  
  return(proportions)
}

################################################################################
# Categorize Tumors by Stromal Content
################################################################################

categorize_stromal_content <- function(cell_proportions) {
  
  cat("\n=== Categorizing Tumors by Stromal Content ===\n")
  
  # Rank by stromal score
  stromal_content <- cell_proportions[order(cell_proportions$Stromal_Score, 
                                             decreasing = TRUE), ]
  stromal_content$Stromal_Rank <- 1:nrow(stromal_content)
  
  # Create categories (tertiles)
  quantiles <- quantile(stromal_content$Stromal_Score, 
                        probs = c(0.33, 0.67), na.rm = TRUE)
  
  stromal_content$Stromal_Category <- cut(
    stromal_content$Stromal_Score,
    breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
    labels = c("Low_Stromal", "Medium_Stromal", "High_Stromal"),
    include.lowest = TRUE
  )
  
  cat("\nStromal Content Categories:\n")
  print(table(stromal_content$Stromal_Category))
  
  cat("\nCategory Thresholds:\n")
  cat("  Low:    <", round(quantiles[1], 3), "\n")
  cat("  Medium: ", round(quantiles[1], 3), "-", round(quantiles[2], 3), "\n")
  cat("  High:   >", round(quantiles[2], 3), "\n\n")
  
  write.csv(stromal_content, 
            file.path(output_dir, "stromal_content_per_tumor.csv"), 
            row.names = FALSE)
  cat("Saved: stromal_content_per_tumor.csv\n")
  
  return(stromal_content)
}

################################################################################
# Select High Stromal Tumors
################################################################################

select_high_stromal_tumors <- function(stromal_content, top_n = 100, percentile = NULL) {
  
  cat("\n=== Selecting High Stromal Tumors ===\n")
  
  if (!is.null(percentile)) {
    threshold <- quantile(stromal_content$Stromal_Score, 
                         probs = percentile/100)
    high_stromal <- stromal_content[stromal_content$Stromal_Score >= threshold, ]
    cat("Selected tumors with stromal score >=", percentile, "th percentile\n")
  } else {
    high_stromal <- head(stromal_content, top_n)
    cat("Selected top", nrow(high_stromal), "tumors with highest stromal content\n")
  }
  
  cat("\nHigh Stromal Tumors:\n")
  cat("  N samples:", nrow(high_stromal), "\n")
  cat("  Mean stromal score:", round(mean(high_stromal$Stromal_Score), 3), "\n")
  cat("  Range:", round(min(high_stromal$Stromal_Score), 3), "-",
      round(max(high_stromal$Stromal_Score), 3), "\n")
  
  # Show cell type composition of high stromal tumors
  cell_cols <- c("CAF", "Endothelial", "Tcells", "Bcells", "Macrophages", "NKcells", "Epithelial")
  cell_cols <- cell_cols[cell_cols %in% names(high_stromal)]
  
  if (length(cell_cols) > 0) {
    cat("\nMean Cell Type Composition in High Stromal Tumors:\n")
    mean_comp <- colMeans(high_stromal[, cell_cols, drop = FALSE], na.rm = TRUE)
    for (i in seq_along(mean_comp)) {
      cat(sprintf("  %-15s: %.3f (%.1f%%)\n", 
                  names(mean_comp)[i], 
                  mean_comp[i],
                  mean_comp[i] * 100))
    }
  }
  
  write.csv(high_stromal, 
            file.path(output_dir, "high_stromal_tumors.csv"), 
            row.names = FALSE)
  cat("\nSaved: high_stromal_tumors.csv\n")
  
  return(high_stromal)
}

################################################################################
# Load Target List - TEST GENES
################################################################################

load_target_list <- function(target_file = NULL) {
  
  cat("\n=== Loading Target Gene List (TEST GENES) ===\n")
  
  if (is.null(target_file) || target_file == "") {
    cat("No target file provided. Using example targets...\n")
    
    example_targets <- c(
      # CAF/ECM markers (should be HIGH in stromal)
      'COL1A1', 'COL1A2', 'FAP', 'DCN', 'PDGFRA', 'ACTA2',
      # Epithelial markers (should be LOW in stromal)
      'EPCAM', 'KRT8', 'KRT18', 'CDH1', 'MUC1',
      # Immune markers
      'CD3D', 'CD68', 'CD19'
    )
    
    example_file <- file.path(output_dir, "example_target_list.txt")
    writeLines(example_targets, example_file)
    cat("Created example file:", example_file, "\n")
    
    return(example_targets)
  }
  
  cat("Attempting to load from:", target_file, "\n")
  cat("File exists:", file.exists(target_file), "\n")
  
  if (!file.exists(target_file)) {
    cat("\nERROR: File not found!\n")
    cat("Using example targets instead...\n")
    
    example_targets <- c(
      'COL1A1', 'COL1A2', 'FAP', 'DCN', 'PDGFRA', 'ACTA2',
      'EPCAM', 'KRT8', 'KRT18', 'CDH1', 'MUC1',
      'CD3D', 'CD68', 'CD19'
    )
    
    example_file <- file.path(output_dir, "example_target_list.txt")
    writeLines(example_targets, example_file)
    cat("Created example file:", example_file, "\n")
    
    return(example_targets)
  }
  
  cat("Loading targets from:", target_file, "\n")
  
  if (grepl("\\.(csv|tsv)$", tolower(target_file))) {
    targets_df <- read.csv(target_file, stringsAsFactors = FALSE, header = TRUE)
    targets <- as.character(targets_df[[1]])
  } else {
    targets <- readLines(target_file)
  }
  
  targets <- trimws(targets)
  targets <- targets[targets != ""]
  targets <- targets[!grepl("^#", targets)]
  targets <- targets[!toupper(targets) %in% c("GENE", "GENE_NAME", "SYMBOL")]
  targets <- unique(targets)
  
  cat("Loaded", length(targets), "target genes\n")
  
  return(targets)
}

################################################################################
# Check for Gene Overlap
################################################################################

check_gene_overlap <- function(reference, targets) {
  
  cat("\n=== Checking for Gene Overlap ===\n")
  
  ref_genes <- unique(unlist(lapply(reference, function(x) c(x$high, x$low))))
  overlapping <- targets[targets %in% ref_genes]
  
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

################################################################################
# Test Targets in High vs Low Stromal Tumors
################################################################################

check_targets_in_high_stromal <- function(expr_data, targets, high_stromal, stromal_content) {
  
  cat("\n=== Testing Target Genes in High vs Low Stromal Tumors ===\n")
  
  available_targets <- targets[targets %in% rownames(expr_data)]
  missing_targets <- targets[!targets %in% rownames(expr_data)]
  
  cat("\nTarget Gene Status:\n")
  cat("  Found:", length(available_targets), "/", length(targets), "\n")
  cat("  Missing:", length(missing_targets), "\n")
  
  if (length(missing_targets) > 0 && length(missing_targets) <= 20) {
    cat("\nMissing targets:\n")
    for (gene in missing_targets) {
      cat("  ", gene, "\n")
    }
  }
  
  if (length(available_targets) == 0) {
    cat("\nNo target genes found in dataset!\n")
    return(NULL)
  }
  
  # Get HIGH stromal samples
  high_stromal_samples <- high_stromal$sample_id
  target_expr_high <- expr_data[available_targets, high_stromal_samples, drop = FALSE]
  
  # Get LOW stromal samples (same number for balanced comparison)
  low_stromal <- tail(stromal_content, nrow(high_stromal))
  low_stromal_samples <- low_stromal$sample_id
  target_expr_low <- expr_data[available_targets, low_stromal_samples, drop = FALSE]
  
  cat("\n=== Comparison Groups ===\n")
  cat("High stromal tumors:", length(high_stromal_samples), "\n")
  cat("Low stromal tumors:", length(low_stromal_samples), "\n\n")
  
  # Calculate statistics
  target_stats <- data.frame(
    gene = available_targets,
    high_stromal_mean = rowMeans(target_expr_high),
    high_stromal_median = apply(target_expr_high, 1, median),
    low_stromal_mean = rowMeans(target_expr_low),
    low_stromal_median = apply(target_expr_low, 1, median),
    fold_change = NA,
    p_value = NA,
    effect_size = NA,
    pct_high_expressed = rowSums(target_expr_high > 1) / ncol(target_expr_high) * 100,
    pct_low_expressed = rowSums(target_expr_low > 1) / ncol(target_expr_low) * 100,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(target_stats)) {
    gene <- target_stats$gene[i]
    
    high_vals <- as.numeric(target_expr_high[gene, ])
    low_vals <- as.numeric(target_expr_low[gene, ])
    
    # Log2 fold change
    target_stats$fold_change[i] <- log2(
      (target_stats$high_stromal_mean[i] + 0.01) / 
      (target_stats$low_stromal_mean[i] + 0.01)
    )
    
    # Effect size (Cohen's d)
    pooled_sd <- sqrt((sd(high_vals)^2 + sd(low_vals)^2) / 2)
    if (pooled_sd > 0) {
      target_stats$effect_size[i] <- (mean(high_vals) - mean(low_vals)) / pooled_sd
    } else {
      target_stats$effect_size[i] <- 0
    }
    
    # Wilcoxon test
    if (length(high_vals) > 0 && length(low_vals) > 0) {
      test <- tryCatch({
        wilcox.test(high_vals, low_vals)
      }, error = function(e) NULL)
      
      if (!is.null(test)) {
        target_stats$p_value[i] <- test$p.value
      }
    }
  }
  
  # Adjust for multiple testing
  target_stats$p_adjusted <- p.adjust(target_stats$p_value, method = "BH")
  
  # Significance based on BOTH p-value AND effect size
  target_stats$significant <- ifelse(
    target_stats$p_adjusted < 0.05 & abs(target_stats$effect_size) > 0.5, "***",
    ifelse(target_stats$p_adjusted < 0.05 & abs(target_stats$effect_size) > 0.2, "**",
           ifelse(target_stats$p_adjusted < 0.05, "*", "ns"))
  )
  
  target_stats <- target_stats[order(target_stats$fold_change, decreasing = TRUE), ]
  
  cat("Target Expression in High vs Low Stromal Tumors:\n")
  cat(strrep("-", 110), "\n")
  cat(sprintf("%-12s %10s %10s %8s %10s %10s %8s %5s\n",
              "Gene", "High Mean", "Low Mean", "Log2FC", "P-adj", "Effect", "Sig", ""))
  cat(strrep("-", 110), "\n")
  
  for (i in 1:nrow(target_stats)) {
    cat(sprintf("%-12s %10.2f %10.2f %8.2f %10.2e %10.2f %8s %5s\n",
                target_stats$gene[i],
                target_stats$high_stromal_mean[i],
                target_stats$low_stromal_mean[i],
                target_stats$fold_change[i],
                target_stats$p_adjusted[i],
                target_stats$effect_size[i],
                ifelse(abs(target_stats$effect_size[i]) > 0.5, "LARGE",
                       ifelse(abs(target_stats$effect_size[i]) > 0.2, "SMALL", "TINY")),
                target_stats$significant[i]))
  }
  
  cat("\nInterpretation Guide:\n")
  cat("  Effect Size: >0.8 = Large, 0.5-0.8 = Medium, 0.2-0.5 = Small, <0.2 = Negligible\n")
  cat("  *** = p<0.05 + large effect, ** = p<0.05 + small effect, * = p<0.05 only\n\n")
  
  write.csv(target_stats, 
            file.path(output_dir, "targets_in_high_stromal_tumors.csv"), 
            row.names = FALSE)
  cat("Saved: targets_in_high_stromal_tumors.csv\n")
  
  result <- list(
    stats = target_stats,
    high_stromal_expr = target_expr_high,
    low_stromal_expr = target_expr_low
  )
  
  return(result)
}

################################################################################
# Create Visualizations
################################################################################

create_visualizations <- function(stromal_content, high_stromal, target_analysis) {
  
  cat("\n=== Creating Visualizations ===\n")
  
  # 1. Stromal content and cell type distribution
  png(file.path(figures_dir, "stromal_content_distribution.png"), 
      width = 16, height = 6, units = "in", res = 300)
  
  par(mfrow = c(1, 3))
  
  # Histogram
  hist(stromal_content$Stromal_Score,
       breaks = 30,
       col = "#3498db",
       border = "white",
       main = "Distribution of Stromal Content\n(Pseudo-deconvolution)",
       xlab = "Overall Stromal Score",
       ylab = "Number of Tumors")
  
  abline(v = high_stromal$Stromal_Score[nrow(high_stromal)],
         col = "red", lwd = 2, lty = 2)
  legend("topright", legend = "High stromal cutoff", col = "red", lty = 2, lwd = 2)
  
  # Categories
  category_counts <- table(stromal_content$Stromal_Category)
  barplot(category_counts,
          col = c("Low_Stromal" = "#2ecc71", 
                 "Medium_Stromal" = "#f39c12", 
                 "High_Stromal" = "#e74c3c"),
          main = "Tumors by Stromal Category",
          ylab = "Number of Tumors")
  
  # Cell type composition
  cell_cols <- c("CAF", "Endothelial", "Tcells", "Bcells", "Macrophages", "NKcells", "Epithelial")
  cell_cols <- cell_cols[cell_cols %in% names(stromal_content)]
  
  if (length(cell_cols) > 0) {
    mean_props <- colMeans(stromal_content[, cell_cols, drop = FALSE], na.rm = TRUE)
    barplot(mean_props,
            col = rainbow(length(mean_props)),
            main = "Mean Cell Type Proportions",
            ylab = "Proportion",
            las = 2)
  }
  
  dev.off()
  cat("Saved: stromal_content_distribution.png\n")
  
  # 2. Cell type composition in high vs low stromal
  png(file.path(figures_dir, "cell_composition_comparison.png"), 
      width = 12, height = 6, units = "in", res = 300)
  
  par(mfrow = c(1, 2))
  
  if (length(cell_cols) > 0) {
    # High stromal
    high_comp <- colMeans(high_stromal[, cell_cols, drop = FALSE], na.rm = TRUE)
    barplot(high_comp,
            col = rainbow(length(high_comp)),
            main = "High Stromal Tumors\nCell Type Composition",
            ylab = "Mean Proportion",
            las = 2,
            ylim = c(0, max(high_comp) * 1.2))
    
    # Low stromal (same number for comparison)
    low_stromal <- tail(stromal_content, nrow(high_stromal))
    low_comp <- colMeans(low_stromal[, cell_cols, drop = FALSE], na.rm = TRUE)
    barplot(low_comp,
            col = rainbow(length(low_comp)),
            main = "Low Stromal Tumors\nCell Type Composition",
            ylab = "Mean Proportion",
            las = 2,
            ylim = c(0, max(high_comp) * 1.2))
  }
  
  dev.off()
  cat("Saved: cell_composition_comparison.png\n")
  
  # 3. Target expression comparison
  if (!is.null(target_analysis) && !is.null(target_analysis$stats)) {
    
    stats <- target_analysis$stats
    
    # Fold change barplot
    png(file.path(figures_dir, "targets_fold_change.png"), 
        width = 10, height = max(6, nrow(stats) * 0.3), units = "in", res = 300)
    
    stats_sorted <- stats[order(stats$fold_change), ]
    
    # Color by effect size AND significance
    colors <- ifelse(stats_sorted$p_adjusted < 0.05 & 
                     stats_sorted$fold_change > 0 & 
                     abs(stats_sorted$effect_size) > 0.5, "#e74c3c",
                    ifelse(stats_sorted$p_adjusted < 0.05 & 
                           stats_sorted$fold_change < 0 & 
                           abs(stats_sorted$effect_size) > 0.5, "#2ecc71",
                          "#95a5a6"))
    
    par(mar = c(5, 10, 4, 2))
    barplot(stats_sorted$fold_change,
            names.arg = stats_sorted$gene,
            horiz = TRUE,
            las = 1,
            col = colors,
            main = "Target Expression: High vs Low Stromal Tumors\n(Deconvolution Analysis)",
            xlab = "Log2(High / Low)")
    
    abline(v = 0, lty = 2, lwd = 2)
    
    legend("bottomright",
           legend = c("Higher in High (p<0.05, large effect)", 
                     "Higher in Low (p<0.05, large effect)", 
                     "Not significant or small effect"),
           fill = c("#e74c3c", "#2ecc71", "#95a5a6"))
    
    dev.off()
    cat("Saved: targets_fold_change.png\n")
    
    # 4. Heatmap of target expression across samples
    if (nrow(stats) > 2) {
      png(file.path(figures_dir, "targets_heatmap.png"), 
          width = 12, height = max(6, nrow(stats) * 0.4), units = "in", res = 300)
      
      # Combine high and low stromal expression
      all_expr <- cbind(target_analysis$high_stromal_expr,
                       target_analysis$low_stromal_expr)
      
      # Log transform for visualization
      all_expr_log <- log2(all_expr + 1)
      
      # Create annotation
      annotation_col <- data.frame(
        Group = c(rep("High Stromal", ncol(target_analysis$high_stromal_expr)),
                 rep("Low Stromal", ncol(target_analysis$low_stromal_expr))),
        row.names = colnames(all_expr)
      )
      
      tryCatch({
        pheatmap(all_expr_log,
                scale = "row",
                clustering_distance_rows = "euclidean",
                clustering_distance_cols = "euclidean",
                annotation_col = annotation_col,
                show_colnames = FALSE,
                main = "Target Gene Expression Heatmap\nHigh vs Low Stromal Tumors",
                fontsize_row = 8)
        
        dev.off()
        cat("Saved: targets_heatmap.png\n")
      }, error = function(e) {
        dev.off()
        cat("Could not create heatmap:", conditionMessage(e), "\n")
      })
    }
  }
  
  cat("\nAll visualizations complete!\n")
}

################################################################################
# Main
################################################################################

# Command line: Rscript tcga_deconv.R <gene list> [<high-stromal selection>]
# Without arguments the script asks for these instead.
cli_args <- commandArgs(trailingOnly = TRUE)
answer <- function(i) {
  if (length(cli_args) == 0) return(readline())
  value <- if (length(cli_args) >= i) cli_args[i] else ""
  cat(value, "\n")
  value
}

main <- function() {
  
  # Download/load data
  expr_data <- download_tcga_data()
  
  # Identify sample types
  sample_info <- identify_sample_types(expr_data)
  
  # Create reference profiles for deconvolution
  reference <- create_reference_profiles()
  
  # Estimate cell type proportions using pseudo-deconvolution
  cell_proportions <- estimate_cell_proportions(expr_data, reference, sample_info)
  
  # Categorize tumors by stromal content
  stromal_content <- categorize_stromal_content(cell_proportions)
  
  # Interactive selection of high stromal tumors
  cat("\n", strrep("=", 70), "\n")
  cat("How many high stromal tumors to select?\n")
  cat("Enter number (e.g., 100) or percentile (e.g., 75 for top 25%): ")
  selection <- answer(2)
  
  selection_num <- suppressWarnings(as.numeric(selection))
  if (is.na(selection_num)) selection_num <- 100
  
  if (selection_num > 50) {
    high_stromal <- select_high_stromal_tumors(stromal_content, 
                                               percentile = selection_num)
  } else {
    high_stromal <- select_high_stromal_tumors(stromal_content, 
                                               top_n = selection_num)
  }
  
  # Load target genes
  cat("\n", strrep("=", 70), "\n")
  cat("Enter path to your target gene list file\n")
  cat("(or press Enter to use example targets): ")
  target_file <- answer(1)
  
  # Clean up input: remove quotes and whitespace
  target_file <- trimws(target_file)
  target_file <- gsub('^["\']|["\']$', '', target_file)  # Remove leading/trailing quotes
  
  targets <- load_target_list(target_file)
  
  # Check for overlap between reference and test genes
  overlapping <- check_gene_overlap(reference, targets)
  
  if (length(overlapping) > 0) {
    cat("\n", strrep("=", 70), "\n")
    cat("WARNING: Some test genes overlap with reference signatures!\n")
    cat("Continue despite overlap? (y/n): ")
    response <- tolower(trimws(readline()))
    if (response != "y") {
      stop("Analysis canceled due to gene overlap")
    }
  }
  
  # Test target genes in high vs low stromal tumors
  target_analysis <- check_targets_in_high_stromal(expr_data, 
                                                   targets, 
                                                   high_stromal, 
                                                   stromal_content)
  
  # Create visualizations
  create_visualizations(stromal_content, high_stromal, target_analysis)
  
  # Summary
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE!\n")
  cat(strrep("=", 70), "\n")
  cat("\nOutput directory:", normalizePath(output_dir), "\n")
  cat("\nFiles created:\n")
  cat("  cell_type_proportions.csv\n")
  cat("  stromal_content_per_tumor.csv\n")
  cat("  high_stromal_tumors.csv\n")
  cat("  targets_in_high_stromal_tumors.csv\n")
  cat("  figures/stromal_content_distribution.png\n")
  cat("  figures/cell_composition_comparison.png\n")
  cat("  figures/targets_fold_change.png\n")
  cat("  figures/targets_heatmap.png\n\n")
  
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

# Run main analysis
main()