#!/usr/bin/env R
# Metabric-BRCA stromal component analysis: deconvolution and non-circular enrichment
# ------------------------------------------------------------------------------
# Design:
#   1) The target gene list is removed from the cell-type marker lists used to
#      estimate cell proportions and Stromal_Score (avoids circularity in the score).
#   2) The enrichment analysis excludes all reference marker genes:
#        - reference marker genes are excluded from the ranked background
#        - enrichment is evaluated only on the remaining "novel" ranked genes
#      (Targets are not filtered against the reference markers again; after step 1
#       they should not overlap, and any remaining overlap is reported.)
#   3) Several enrichment plots are drawn in addition to the running-sum curve.
#
# Outputs:
#   - enrichment_ranked_genes.csv
#   - target_enrichment_gsea.csv
#   - figures/target_enrichment_curve.png
#   - figures/target_barcode.png
#   - figures/target_ecdf.png
#   - figures/target_rho_distribution.png
#   - figures/top_targets_lollipop.png
#   - figures/null_es_distribution.png

packages <- c("dplyr", "ggplot2", "pheatmap", "reshape2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

set.seed(1)

################################################################################
# Setup
################################################################################

output_dir <- "./Metabric_stromal_analysis"
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("Metabric-BRCA STROMAL CONTENT ANALYZER\n")
cat("DECONVOLUTION + NON-CIRCULAR ENRICHMENT (v2)\n")
cat("========================================\n\n")
cat("Output directory:", normalizePath(output_dir), "\n")

################################################################################
# Utility helpers
################################################################################

sanitize_gene_symbols <- function(x) {
  x <- trimws(x)
  x <- x[x != ""]
  x <- x[!grepl("^#", x)]
  x <- gsub("\\s+", "", x)
  toupper(unique(x))
}

get_reference_genes <- function(reference) {
  unique(toupper(unlist(lapply(reference, function(ct) c(ct$high, ct$low)))))
}

remove_genes_from_reference <- function(reference, genes_to_remove) {
  genes_to_remove <- toupper(genes_to_remove)

  removed_counts <- data.frame(
    cell_type = character(0),
    removed_high = integer(0),
    removed_low = integer(0),
    stringsAsFactors = FALSE
  )

  for (ct in names(reference)) {
    high0 <- toupper(reference[[ct]]$high)
    low0  <- toupper(reference[[ct]]$low)

    high1 <- setdiff(high0, genes_to_remove)
    low1  <- setdiff(low0, genes_to_remove)

    removed_counts <- rbind(
      removed_counts,
      data.frame(
        cell_type = ct,
        removed_high = length(setdiff(high0, high1)),
        removed_low  = length(setdiff(low0, low1)),
        stringsAsFactors = FALSE
      )
    )

    reference[[ct]]$high <- high1
    reference[[ct]]$low  <- low1
  }

  list(reference = reference, removed = removed_counts)
}

#############################
# Load Metabric expression/metadata
#############################
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

################################################################################
# Load Target List
################################################################################

load_target_list <- function(target_file = NULL) {

  cat("\n=== Loading Target Gene List (TEST GENES) ===\n")

  cat("Attempting to load from:", target_file, "\n")
  cat("File exists:", file.exists(target_file), "\n")

  if (grepl("\\.(csv|tsv)$", tolower(target_file))) {
    targets_df <- read.csv(target_file, stringsAsFactors = FALSE, header = TRUE)
    targets <- as.character(targets_df[[1]])
  } else {
    targets <- readLines(target_file)
  }

  targets <- sanitize_gene_symbols(targets)
  targets <- targets[!targets %in% c("GENE", "GENE_NAME", "SYMBOL")]
  print(targets)
  cat("Loaded", length(targets), "target genes\n")
  return(targets)
}

################################################################################
# Create Reference Cell Type Profiles (Pseudo-deconvolution)
################################################################################

create_reference_profiles <- function() {

  cat("\n=== Creating Cell Type Reference Profiles ===\n")
  cat("Using pseudo-deconvolution with positive/negative markers\n\n")

  reference <- list(

    CAF = list(
      high = c('COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A1', 'COL6A2',
               'FAP', 'ACTA2', 'PDGFRA', 'PDGFRB', 'LUM', 'DCN', 'POSTN',
               'THBS2', 'FN1', 'VCAN', 'BGN', 'SPARC', 'MMP2', 'MMP11'),
      low = c('PTPRC', 'CD3D', 'CD19', 'EPCAM', 'KRT8', 'KRT18', 'CD34', 'PECAM1'),
      weight = 2.0
    ),

    Endothelial = list(
      high = c('PECAM1', 'CD34', 'VWF', 'CDH5', 'ENG', 'PLVAP', 'CLDN5',
               'TEK', 'FLT1', 'KDR', 'ESAM', 'MMRN1', 'CLEC14A'),
      low = c('COL1A1', 'EPCAM', 'CD3D', 'FAP', 'ACTA2', 'CD68'),
      weight = 1.5
    ),

    Tcells = list(
      high = c('CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'IL7R',
               'CD2', 'CD5', 'CD7', 'TRAC', 'TRBC1', 'TRBC2'),
      low = c('CD19', 'CD68', 'EPCAM', 'COL1A1', 'CD34'),
      weight = 0.3
    ),

    Bcells = list(
      high = c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'BLK', 'CD22', 'BANK1'),
      low = c('CD3D', 'CD68', 'EPCAM', 'COL1A1', 'CD34'),
      weight = 0.3
    ),

    Macrophages = list(
      high = c('CD68', 'CD163', 'CD14', 'FCGR3A', 'ITGAX', 'CD86',
               'CSF1R', 'MARCO', 'MSR1', 'MRC1'),
      low = c('CD3D', 'CD19', 'EPCAM', 'KRT8', 'COL1A1'),
      weight = 0.4
    ),

    NKcells = list(
      high = c('NCAM1', 'NCR1', 'KLRB1', 'KLRD1', 'NKG7', 'GNLY'),
      low = c('CD3D', 'CD19', 'CD68', 'EPCAM', 'COL1A1'),
      weight = 0.2
    ),

    Epithelial = list(
      high = c('EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'MUC1',
               'KRT7', 'CLDN4', 'CLDN3', 'TJP1'),
      low = c('PTPRC', 'COL1A1', 'VWF', 'CD3D', 'FAP'),
      weight = -1.0
    )
  )

  for (ct in names(reference)) {
    reference[[ct]]$high <- sanitize_gene_symbols(reference[[ct]]$high)
    reference[[ct]]$low  <- sanitize_gene_symbols(reference[[ct]]$low)
  }

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

estimate_cell_proportions <- function(expr_data, reference) {

  cat("\n=== Running Pseudo-Deconvolution ===\n")
  cat("Estimating cell type proportions using positive/negative markers\n\n")

  tumor_samples <- colnames(expr_data)
  tumor_expr <- expr_data[, tumor_samples, drop = FALSE]

  cat("Analyzing", ncol(tumor_expr), "tumor samples\n\n")

  proportions <- data.frame(
    sample_id = tumor_samples,
    stringsAsFactors = FALSE
  )

  for (cell_type in names(reference)) {

    profile <- reference[[cell_type]]

    high_genes <- profile$high[profile$high %in% rownames(tumor_expr)]
    low_genes  <- profile$low[profile$low %in% rownames(tumor_expr)]

    cat(cell_type, ":\n")
    cat("  High markers found:", length(high_genes), "/", length(profile$high), "\n")
    cat("  Low  markers found:", length(low_genes),  "/", length(profile$low),  "\n")

    if (length(high_genes) == 0) {
      proportions[[cell_type]] <- NA
      cat("  WARNING: No high markers found!\n\n")
      next
    }

    high_expr <- colMeans(tumor_expr[high_genes, , drop = FALSE])

    if (length(low_genes) > 0) {
      low_expr <- colMeans(tumor_expr[low_genes, , drop = FALSE])
      score <- high_expr - low_expr
    } else {
      score <- high_expr
    }

    denom <- (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
    if (is.na(denom) || denom == 0) {
      score <- rep(0, length(score))
    } else {
      score <- (score - min(score, na.rm = TRUE)) / denom
    }

    proportions[[cell_type]] <- score

    cat("  Score range:", round(min(score, na.rm = TRUE), 3),
        "to", round(max(score, na.rm = TRUE), 3), "\n\n")
  }

  cell_cols <- intersect(names(reference), names(proportions))
  row_sums <- rowSums(proportions[, cell_cols, drop = FALSE], na.rm = TRUE)
  row_sums[row_sums == 0 | is.na(row_sums)] <- NA
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
  
  cat("=== Calculating Weighted Stromal Score ===\n")
  proportions$Stromal_Score <- 0
  for (cell_type in names(reference)) {
    if (cell_type %in% names(proportions)) {
      weight <- reference[[cell_type]]$weight
      proportions$Stromal_Score <- proportions$Stromal_Score + (proportions[[cell_type]] * weight)
      cat(sprintf("  %s: weight = %.1f\n", cell_type, weight))
    }
  }

  ss <- proportions$Stromal_Score
  denom <- max(ss, na.rm = TRUE) - min(ss, na.rm = TRUE)
  if (!is.na(denom) && denom > 0) {
    proportions$Stromal_Score <- (ss - min(ss, na.rm = TRUE)) / denom
  } else {
    proportions$Stromal_Score <- 0
  }

  cat("\nStromal Score Statistics:\n")
  cat("  Min:", round(min(proportions$Stromal_Score), 3), "\n")
  cat("  Median:", round(median(proportions$Stromal_Score), 3), "\n")
  cat("  Mean:", round(mean(proportions$Stromal_Score), 3), "\n")
  cat("  Max:", round(max(proportions$Stromal_Score), 3), "\n\n")

  write.csv(proportions, file.path(output_dir, "cell_type_proportions.csv"), row.names = FALSE)
  cat("Saved: cell_type_proportions.csv\n")

  return(proportions)
}

################################################################################
# Categorize Tumors by Stromal Content
################################################################################

categorize_stromal_content <- function(cell_proportions) {

  cat("\n=== Categorizing Tumors by Stromal Content ===\n")

  stromal_content <- cell_proportions[order(cell_proportions$Stromal_Score, decreasing = TRUE), ]
  stromal_content$Stromal_Rank <- 1:nrow(stromal_content)

  quantiles <- quantile(stromal_content$Stromal_Score, probs = c(0.33, 0.67), na.rm = TRUE)

  stromal_content$Stromal_Category <- cut(
    stromal_content$Stromal_Score,
    breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
    labels = c("Low_Stromal", "Medium_Stromal", "High_Stromal"),
    include.lowest = TRUE
  )

  cat("\nStromal Content Categories:\n")
  print(table(stromal_content$Stromal_Category))

  write.csv(stromal_content, file.path(output_dir, "stromal_content_per_tumor.csv"), row.names = FALSE)
  cat("Saved: stromal_content_per_tumor.csv\n")

  stromal_content
}

################################################################################
# Select High Stromal Tumors (kept for continuity)
################################################################################

select_high_stromal_tumors <- function(stromal_content, top_n = 100, percentile = NULL) {

  cat("\n=== Selecting High Stromal Tumors ===\n")

  if (!is.null(percentile)) {
    threshold <- quantile(stromal_content$Stromal_Score, probs = percentile/100, na.rm = TRUE)
    high_stromal <- stromal_content[stromal_content$Stromal_Score >= threshold, ]
    cat("Selected tumors with stromal score >=", percentile, "th percentile\n")
  } else {
    high_stromal <- head(stromal_content, top_n)
    cat("Selected top", nrow(high_stromal), "tumors with highest stromal content\n")
  }

  write.csv(high_stromal, file.path(output_dir, "high_stromal_tumors.csv"), row.names = FALSE)
  cat("Saved: high_stromal_tumors.csv\n")

  high_stromal
}

################################################################################
# Rank ALL "novel" genes by effect size vs Stromal_Score
################################################################################

rank_genes_by_stromal_effect <- function(expr_data, stromal_content, exclude_genes = character(0)) {

  tumor_samples <- stromal_content$sample_id
  stromal_score <- stromal_content$Stromal_Score
  names(stromal_score) <- tumor_samples

  tumor_expr <- expr_data[, tumor_samples, drop = FALSE]

  exclude_genes <- unique(toupper(exclude_genes))
  genes <- setdiff(rownames(tumor_expr), exclude_genes)

  cat("\n=== Ranking genes by effect size vs Stromal_Score ===\n")
  cat("Tumor samples:", length(tumor_samples), "\n")
  cat("Genes (expressed):", nrow(tumor_expr), "\n")
  cat("Excluded genes from ranking:", length(exclude_genes), "\n")
  cat("Genes ranked (novel):", length(genes), "\n\n")

  X <- log2(tumor_expr[genes, , drop = FALSE] + 1)

  rho <- apply(X, 1, function(v) suppressWarnings(cor(v, stromal_score, method = "spearman", use = "pairwise.complete.obs")))
  rho[is.na(rho)] <- 0

  ranked <- sort(rho, decreasing = TRUE)

  ranked_df <- data.frame(gene = names(ranked), spearman_rho = as.numeric(ranked), stringsAsFactors = FALSE)
  write.csv(ranked_df, file.path(output_dir, "enrichment_ranked_genes.csv"), row.names = FALSE)
  cat("Saved: enrichment_ranked_genes.csv\n")

  ranked
}

################################################################################
# Lightweight GSEA (running-sum) + permutation p-value
################################################################################

gsea_enrichment <- function(ranked_stats, gene_set, nperm = 5000, weight_p = 1) {

  all_genes <- names(ranked_stats)
  gene_set <- toupper(gene_set)
  gene_set <- intersect(gene_set, all_genes)

  N <- length(all_genes)
  Nh <- length(gene_set)

  if (Nh < 2) {
    return(list(
      ES = NA, NES = NA, p_value = NA, n_in_set = Nh,
      warning = "Gene set too small after filtering."
    ))
  }

  hits <- all_genes %in% gene_set
  stats_abs_p <- abs(ranked_stats)^weight_p

  NR <- sum(stats_abs_p[hits])

  P_hit <- ifelse(hits, stats_abs_p / NR, 0)
  P_miss <- ifelse(!hits, 1 / (N - Nh), 0)
  running <- cumsum(P_hit - P_miss)

  ES <- running[which.max(abs(running))]

  perm_ES <- numeric(nperm)
  for (i in seq_len(nperm)) {
    perm_set <- sample(all_genes, Nh, replace = FALSE)
    hits_p <- all_genes %in% perm_set
    NRp <- sum(stats_abs_p[hits_p])
    P_hit_p <- ifelse(hits_p, stats_abs_p / NRp, 0)
    P_miss_p <- ifelse(!hits_p, 1 / (N - Nh), 0)
    running_p <- cumsum(P_hit_p - P_miss_p)
    perm_ES[i] <- running_p[which.max(abs(running_p))]
  }

  pval <- (sum(abs(perm_ES) >= abs(ES)) + 1) / (nperm + 1)

  same_sign <- perm_ES[sign(perm_ES) == sign(ES)]
  if (length(same_sign) < 10) {
    NES <- ES / mean(abs(perm_ES))
  } else {
    NES <- ES / mean(abs(same_sign))
  }

  list(
    ES = ES,
    NES = NES,
    p_value = pval,
    n_in_set = Nh,
    running = running,
    hits = hits,
    perm_ES = perm_ES
  )
}

################################################################################
# Enrichment visualization helpers
################################################################################

plot_gsea_curve <- function(gsea_res, title = "Target set enrichment") {
  if (is.null(gsea_res$running) || all(is.na(gsea_res$running))) return(invisible(NULL))

  df <- data.frame(
    rank = seq_along(gsea_res$running),
    running_ES = gsea_res$running,
    hit = gsea_res$hits
  )

  ggplot(df, aes(x = rank, y = running_ES)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(data = df[df$hit, ], aes(x = rank, y = running_ES), size = 0.7, alpha = 0.6) +
    labs(
      title = title,
      subtitle = sprintf("ES=%.3f, NES=%.3f, p=%.3g, n=%d",
                         gsea_res$ES, gsea_res$NES, gsea_res$p_value, gsea_res$n_in_set),
      x = "Gene rank (most stromal-aligned → least)",
      y = "Running enrichment score"
    ) +
    theme_minimal()
}

plot_barcode <- function(ranked_stats, gene_set, title = "Target barcode along ranked list") {
  all_genes <- names(ranked_stats)
  gene_set <- intersect(toupper(gene_set), all_genes)

  df <- data.frame(rank = seq_along(all_genes), hit = all_genes %in% gene_set)

  ggplot(df, aes(x = rank, y = 0)) +
    geom_point(data = df[df$hit, ], aes(x = rank, y = 0), size = 0.8, alpha = 0.7) +
    labs(title = title, x = "Gene rank (most stromal-aligned → least)", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

plot_ecdf <- function(ranked_stats, gene_set, title = "ECDF of target ranks") {
  
  all_genes <- names(ranked_stats)
  gene_set  <- intersect(toupper(gene_set), all_genes)
  
  ranks <- match(gene_set, all_genes)
  keep  <- !is.na(ranks)
  
  if (sum(keep) < 2) return(invisible(NULL))
  
  df_targets <- data.frame(
    gene = gene_set[keep],
    rank = ranks[keep],
    stringsAsFactors = FALSE
  )
  
  df_targets <- df_targets[order(df_targets$rank), ]
  df_targets$frac <- seq_len(nrow(df_targets)) / nrow(df_targets)
  
  df_null <- data.frame(
    rank = c(1, length(all_genes)),
    frac = c(0, 1)
  )
  
  ggplot() +
    # Null expectation
    geom_line(
      data = df_null,
      aes(x = rank, y = frac),
      linetype = 2,
      color = "gray50"
    ) +
    # ECDF line
    geom_line(
      data = df_targets,
      aes(x = rank, y = frac),
      linewidth = 1
    ) +
    # Points for each gene
    geom_point(
      data = df_targets,
      aes(x = rank, y = frac),
      size = 2
    ) +
    # Gene labels
    ggrepel::geom_text_repel(
      data = df_targets,
      aes(x = rank, y = frac, label = gene),
      
      size = 3,
      
      # stronger repulsion
      force = 3,
      force_pull = 0.5,
      
      # push text away from points + line
      point.padding = 0.6,
      box.padding = 0.7,
      
      # keep labels off the curve
      min.segment.length = 0,
      segment.alpha = 0.4,
      segment.size = 0.3,
      
      # allow labels to move freely
      direction = "both",
      max.overlaps = Inf
    ) +
    labs(
      title = title,
      subtitle = "Each point = one target gene (labeled by symbol)",
      x = "Gene rank (most stromal-aligned → least)",
      y = "Cumulative fraction of targets"
    ) +
    theme_minimal() + xlim(-2500, 20000)
}

plot_rho_distribution <- function(ranked_stats, gene_set, title = "Effect-size distribution") {
  all_genes <- names(ranked_stats)
  gene_set <- intersect(toupper(gene_set), all_genes)

  df <- data.frame(
    rho = as.numeric(ranked_stats),
    group = ifelse(all_genes %in% gene_set, "Targets", "Other_novel"),
    stringsAsFactors = FALSE
  )

  # Subsample "Other_novel" for plotting if enormous (keeps file size manageable)
  if (sum(df$group == "Other_novel") > 5000) {
    set.seed(1)
    other_idx <- which(df$group == "Other_novel")
    keep_other <- sample(other_idx, 5000)
    df <- rbind(df[df$group == "Targets", ], df[keep_other, ])
  }

  ggplot(df, aes(x = group, y = rho)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    labs(title = title,
         subtitle = "Spearman rho vs Stromal_Score (novel genes only)",
         x = "", y = "Spearman rho") +
    theme_minimal()
}

plot_top_targets_lollipop <- function(ranked_stats, gene_set, top_n = 20, title = "Top target genes by stromal association") {
  all_genes <- names(ranked_stats)
  gene_set <- intersect(toupper(gene_set), all_genes)

  df <- data.frame(gene = gene_set, rho = ranked_stats[gene_set], stringsAsFactors = FALSE)
  df <- df[order(df$rho, decreasing = TRUE), ]
  df <- head(df, top_n)
  df$gene <- factor(df$gene, levels = rev(df$gene))

  ggplot(df, aes(x = gene, y = rho)) +
    geom_segment(aes(x = gene, xend = gene, y = 0, yend = rho)) +
    geom_point(size = 2) +
    coord_flip() +
    labs(title = title, x = "", y = "Spearman rho") +
    theme_minimal()
}

plot_null_es <- function(perm_ES, ES_obs, title = "Null distribution of ES") {
  df <- data.frame(ES = perm_ES)
  ggplot(df, aes(x = ES)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = ES_obs, linetype = 2, linewidth = 1) +
    geom_vline(xintercept = -ES_obs, linetype = 2, linewidth = 1) +
    labs(title = title,
         subtitle = "Dashed lines = ± observed ES",
         x = "Permuted ES", y = "Count") +
    theme_minimal()
}

################################################################################
# Main
################################################################################

# Command line: Rscript metabric_brca_stromal_enrichment.R <gene list> [<high-stromal selection>]
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
  expr_data <- metabric$expression

  cat("\n", strrep("=", 70), "\n")
  cat("Enter path to your target/marker gene list file\n")
  cat("(or press Enter to use example targets): ")
  target_file <- trimws(answer(1))
  target_file <- gsub('^["\']|["\']$', '', target_file)

  targets <- load_target_list(target_file)

  reference <- create_reference_profiles()

  cat("\n=== Removing provided targets from cell-type marker lists (anti-circularity) ===\n")
  removed <- remove_genes_from_reference(reference, targets)
  reference <- removed$reference
  cat("Removed counts by cell type:\n")
  print(removed$removed)
  cat("\n")

  cell_proportions <- estimate_cell_proportions(expr_data, reference) # , sample_info
  stromal_content <- categorize_stromal_content(cell_proportions)

  cat("\n", strrep("=", 70), "\n")
  cat("How many high stromal tumors to select?\n")
  cat("Enter number (e.g., 100) or percentile (e.g., 75 for top 25%): ")
  selection <- answer(2)
  selection_num <- suppressWarnings(as.numeric(selection))
  if (is.na(selection_num)) selection_num <- 100

  if (selection_num > 50) {
    high_stromal <- select_high_stromal_tumors(stromal_content, percentile = selection_num)
  } else {
    high_stromal <- select_high_stromal_tumors(stromal_content, top_n = selection_num)
  }

  cat("\n", strrep("=", 70), "\n")
  cat("=== NON-CIRCULAR ENRICHMENT (GSEA) ===\n")

  # Exclude reference marker genes from the ENTIRE enrichment analysis
  reference_genes_used <- get_reference_genes(reference)

  ranked <- rank_genes_by_stromal_effect(
    expr_data = expr_data,
    stromal_content = stromal_content,
    exclude_genes = reference_genes_used
  )

  # Targets are not filtered against the reference markers again here: they were
  # already removed from the reference marker lists, so the two should not overlap.
  targets_available <- intersect(targets, rownames(expr_data))

  # Report if any overlap still exists (should be rare unless case/alias issues)
  overlap_now <- intersect(targets_available, reference_genes_used)
  if (length(overlap_now) > 0) {
    cat("\nWARNING: Some targets still overlap reference markers even after removal.\n")
    cat("This can happen due to symbol mismatches/aliases. Overlapping genes:\n")
    print(overlap_now)
  } else {
    cat("\nTargets appear non-overlapping with reference marker genes (as expected).\n")
  }

  cat("\nTargets provided:", length(targets), "\n")
  cat("Targets found in expression matrix:", length(targets_available), "\n")
  cat("Reference marker genes excluded from enrichment background:", length(reference_genes_used), "\n")
  cat("Ranked genes (novel background):", length(ranked), "\n\n")

  gsea_res <- gsea_enrichment(ranked_stats = ranked, gene_set = targets_available, nperm = 5000, weight_p = 1)

  enrich_df <- data.frame(
    ES = gsea_res$ES,
    NES = gsea_res$NES,
    p_value = gsea_res$p_value,
    n_targets_in_set = gsea_res$n_in_set,
    n_ranked_genes = length(ranked),
    stringsAsFactors = FALSE
  )

  write.csv(enrich_df, file.path(output_dir, "target_enrichment_gsea.csv"), row.names = FALSE)
  cat("Saved: target_enrichment_gsea.csv\n")
  print(enrich_df)

  # -----------------------------
  # Plots
  # -----------------------------
  # 1) Running-sum curve
  png(file.path(figures_dir, "target_enrichment_curve.png"), width = 11, height = 6, units = "in", res = 300)
  p1 <- plot_gsea_curve(gsea_res, title = "Target enrichment among novel genes ranked by stromal association")
  if (!is.null(p1)) print(p1)
  dev.off()
  cat("Saved: figures/target_enrichment_curve.png\n")

  # 2) Barcode (hit rug)
  png(file.path(figures_dir, "target_barcode.png"), width = 11, height = 2.5, units = "in", res = 300)
  p2 <- plot_barcode(ranked, targets_available, title = "Target barcode along ranked list (novel genes only)")
  if (!is.null(p2)) print(p2)
  dev.off()
  cat("Saved: figures/target_barcode.png\n")

  # 3) ECDF of target ranks
  png(file.path(figures_dir, "target_ecdf.png"), width = 7, height = 5, units = "in", res = 300)
  p3 <- plot_ecdf(ranked, targets_available, title = "ECDF of target ranks (novel genes only)")
  if (!is.null(p3)) print(p3)
  dev.off()
  cat("Saved: figures/target_ecdf.png\n")

  # 4) Effect-size distribution: targets vs other novel genes
  png(file.path(figures_dir, "target_rho_distribution.png"), width = 7, height = 5, units = "in", res = 300)
  p4 <- plot_rho_distribution(ranked, targets_available, title = "Targets vs background effect-size distribution")
  if (!is.null(p4)) print(p4)
  dev.off()
  cat("Saved: figures/target_rho_distribution.png\n")

  # 5) Top target hits lollipop
  png(file.path(figures_dir, "top_targets_lollipop.png"), width = 7, height = 6, units = "in", res = 300)
  p5 <- plot_top_targets_lollipop(ranked, targets_available, top_n = 25, title = "Top target genes by stromal association")
  if (!is.null(p5)) print(p5)
  dev.off()
  cat("Saved: figures/top_targets_lollipop.png\n")

  # 6) Null ES distribution
  if (!is.null(gsea_res$perm_ES)) {
    png(file.path(figures_dir, "null_es_distribution.png"), width = 7, height = 5, units = "in", res = 300)
    p6 <- plot_null_es(gsea_res$perm_ES, gsea_res$ES, title = "Null ES distribution (permutations)")
    if (!is.null(p6)) print(p6)
    dev.off()
    cat("Saved: figures/null_es_distribution.png\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE!\n")
  cat(strrep("=", 70), "\n")
  cat("\nFiles created/updated:\n")
  cat("  cell_type_proportions.csv\n")
  cat("  stromal_content_per_tumor.csv\n")
  cat("  high_stromal_tumors.csv\n")
  cat("  enrichment_ranked_genes.csv\n")
  cat("  target_enrichment_gsea.csv\n")
  cat("  figures/target_enrichment_curve.png\n")
  cat("  figures/target_barcode.png\n")
  cat("  figures/target_ecdf.png\n")
  cat("  figures/target_rho_distribution.png\n")
  cat("  figures/top_targets_lollipop.png\n")
  cat("  figures/null_es_distribution.png\n\n")

  cat("NON-CIRCULARITY:\n")
  cat("  - The target list is removed from the marker signatures used to build Stromal_Score.\n")
  cat("  - All reference marker genes are excluded from the enrichment background ranking.\n")
  cat("  - Enrichment is evaluated only among novel genes (non-marker genes).\n\n")
}

main()
