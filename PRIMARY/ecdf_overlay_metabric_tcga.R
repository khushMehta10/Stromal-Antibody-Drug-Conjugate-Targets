#!/usr/bin/env Rscript
# Overlay ECDFs (same axis) for METABRIC vs TCGA using the same target gene list.
# - Labels target genes ONLY for METABRIC (TCGA unlabeled)
# - Uses neutral (non-black) line colors and different linetypes
#
# Expected inputs:
#   METABRIC: ../Metabric/Metabric_stromal_analysis/enrichment_ranked_genes.csv
#   TCGA:     ../TCGA_Data/tcga_stromal_analysis/enrichment_ranked_genes.csv
#
# The enrichment_ranked_genes.csv is written by:
#   metabric_brca_stromal_enrichment.R and tcga_brca_stromal_enrichment.R
#
# Output:
#   ./combined_ecdf/figures/target_ecdf_overlay_metabric_tcga.png

packages <- c("dplyr", "ggplot2", "ggrepel")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

sanitize_gene_symbols <- function(x) {
  x <- trimws(x)
  x <- x[x != ""]
  x <- x[!grepl("^#", x)]
  x <- gsub("\\s+", "", x)
  toupper(unique(x))
}

load_target_list <- function(target_file) {
  if (is.null(target_file) || target_file == "" || !file.exists(target_file)) {
    stop("Target file not found. Provide a valid path to your target gene list.")
  }

  if (grepl("\\.(csv|tsv)$", tolower(target_file))) {
    targets_df <- read.csv(target_file, stringsAsFactors = FALSE, header = TRUE)
    targets <- as.character(targets_df[[1]])
  } else {
    targets <- readLines(target_file)
  }

  targets <- sanitize_gene_symbols(targets)
  targets <- targets[!targets %in% c("GENE", "GENE_NAME", "SYMBOL")]
  targets
}

# Read ranked stats from the CSV written by the enrichment scripts
# CSV format: columns gene, spearman_rho
read_ranked_from_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Ranked genes file not found at: ", path,
         "\nRun the corresponding pipeline first so it produces enrichment_ranked_genes.csv")
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("gene", "spearman_rho") %in% colnames(df))) {
    stop("Unexpected CSV format at: ", path, "\nExpected columns: gene, spearman_rho")
  }
  df$gene <- toupper(df$gene)
  stats <- df$spearman_rho
  names(stats) <- df$gene
  stats
}

make_ecdf_df <- function(ranked_stats, gene_set, dataset) {
  all_genes <- names(ranked_stats)
  gene_set  <- intersect(toupper(gene_set), all_genes)

  ranks <- match(gene_set, all_genes)
  keep  <- !is.na(ranks)

  if (sum(keep) < 2) return(NULL)

  df <- data.frame(
    gene = gene_set[keep],
    rank = ranks[keep],
    stringsAsFactors = FALSE
  )
  df <- df[order(df$rank), ]
  df$frac <- seq_len(nrow(df)) / nrow(df)
  df$dataset <- dataset
  df
}

plot_ecdf_overlay_metabric_tcga <- function(
  ranked_metabric,
  ranked_tcga,
  gene_set,
  title = "ECDF of target ranks (novel genes only): METABRIC vs TCGA",
  xlim = c(-2500, 20000),
  colors = c("METABRIC" = "slategray4", "TCGA" = "#28bf50") # neutral, not black
) {
  df_m <- make_ecdf_df(ranked_metabric, gene_set, "METABRIC")
  df_t <- make_ecdf_df(ranked_tcga,     gene_set, "TCGA")

  if (is.null(df_m) && is.null(df_t)) return(invisible(NULL))
  df_all <- rbind(df_m, df_t)

  # METABRIC-only labels
  df_labels <- subset(df_all, dataset == "METABRIC")

  # Null expectation diagonal: (1,0) to (N,1) for each dataset's N is different,
  # but we want a single reference line; use the larger N for the x-axis.
  Nmax <- max(length(ranked_metabric), length(ranked_tcga))
  # Stop the line at the right edge of the x-axis (keeping its true slope 1/Nmax);
  # otherwise xlim() drops the end point and the line is not drawn at all.
  x_end <- min(Nmax, xlim[2])
  df_null <- data.frame(rank = c(1, x_end), frac = c(0, x_end / Nmax))

  ggplot(df_all, aes(x = rank, y = frac, group = dataset)) +
    # null expectation
    geom_line(
      data = df_null,
      aes(x = rank, y = frac),
      inherit.aes = FALSE,
      linetype = 2,
      linewidth = 0.7,
      color = "gray60"
    ) +
    # draw TCGA first so METABRIC can sit on top
    geom_line(
      data = subset(df_all, dataset == "TCGA"),
      aes(color = dataset, linetype = dataset),
      linewidth = 1
    ) +
    geom_point(
      data = subset(df_all, dataset == "TCGA"),
      aes(color = dataset),
      size = 1.3,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_line(
      data = subset(df_all, dataset == "METABRIC"),
      aes(color = dataset, linetype = dataset),
      linewidth = 1.1
    ) +
    geom_line(
      data = df_null,
      aes(x = rank, y = frac),
      inherit.aes = FALSE,
      linetype = 2,
      linewidth = 0.7,
      color = "gray60"
    ) +
    geom_point(
      data = df_labels,
      aes(color = dataset),
      size = 2,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = df_labels,
      aes(label = gene),
      color = "black",
      size = 3,
      force = 3,
      force_pull = 0.5,
      point.padding = 0.6,
      box.padding = 0.7,
      min.segment.length = 0,
      segment.alpha = 0.35,
      segment.size = 0.3,
      direction = "both",
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = c("TCGA" = "solid", "METABRIC" = "solid")) +
    labs(
      title = title,
      subtitle = "METABRIC genes labeled; TCGA unlabeled; dashed gray = null expectation",
      x = "Gene rank (most stromal-aligned → least)",
      y = "Cumulative fraction of targets",
      color = "",
      linetype = ""
    ) +
    theme_minimal() +
    coord_cartesian(xlim = xlim) + coord_fixed(ratio = 20000)
}

# Command line: Rscript ecdf_overlay_metabric_tcga.R <gene list>
# Without arguments the script asks for these instead.
cli_args <- commandArgs(trailingOnly = TRUE)
answer <- function(i) {
  if (length(cli_args) == 0) return(readline())
  value <- if (length(cli_args) >= i) cli_args[i] else ""
  cat(value, "\n")
  value
}

main <- function() {
  # Results written by the METABRIC and TCGA enrichment scripts in their own folders:
  metabric_dir <- "../Metabric/Metabric_stromal_analysis"
  tcga_dir     <- "../TCGA_Data/tcga_stromal_analysis"

  metabric_ranked_csv <- file.path(metabric_dir, "enrichment_ranked_genes.csv")
  tcga_ranked_csv     <- file.path(tcga_dir, "enrichment_ranked_genes.csv")

  out_dir <- "./combined_ecdf"
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n=== ECDF OVERLAY: METABRIC vs TCGA ===\n")
  cat("METABRIC ranked CSV:", normalizePath(metabric_ranked_csv, winslash = "/", mustWork = FALSE), "\n")
  cat("TCGA ranked CSV:    ", normalizePath(tcga_ranked_csv, winslash = "/", mustWork = FALSE), "\n\n")

  cat("Enter path to your target gene list file: ")
  target_file <- trimws(answer(1))
  target_file <- gsub('^["\']|["\']$', '', target_file)

  targets <- load_target_list(target_file)

  ranked_metabric <- read_ranked_from_csv(metabric_ranked_csv)
  ranked_tcga     <- read_ranked_from_csv(tcga_ranked_csv)

  # Use only targets present in both ranked backgrounds (both cohorts use the same gene list)
  targets_shared <- intersect(targets, intersect(names(ranked_metabric), names(ranked_tcga)))

  cat("\nTargets provided:              ", length(targets), "\n")
  cat("Targets present in METABRIC:  ", sum(targets %in% names(ranked_metabric)), "\n")
  cat("Targets present in TCGA:      ", sum(targets %in% names(ranked_tcga)), "\n")
  cat("Targets used (shared):        ", length(targets_shared), "\n\n")

  if (length(targets_shared) < 2) {
    stop("Fewer than 2 shared targets were found in both datasets' ranked lists. Check gene symbols/case.")
  }

  out_png <- file.path(fig_dir, "target_ecdf_overlay_metabric_tcga.png")
  png(out_png, width = 7, height = 5, units = "in", res = 300)
  p <- plot_ecdf_overlay_metabric_tcga(
    ranked_metabric = ranked_metabric,
    ranked_tcga     = ranked_tcga,
    gene_set        = targets_shared,
    title           = "ECDF of target ranks (novel genes only): METABRIC vs TCGA",
    xlim            = c(-2500, 20000),
    colors          = c("METABRIC" = "slategray4", "TCGA" = "#28bf50")
  ) + xlim(-2500, 20000)
  if (!is.null(p)) print(p)
  dev.off()

  cat("Saved overlay ECDF to:\n  ", normalizePath(out_png, winslash = "/", mustWork = FALSE), "\n\n")
}

# Run if executed as a script
if (sys.nframe() == 0) {
  main()
}
