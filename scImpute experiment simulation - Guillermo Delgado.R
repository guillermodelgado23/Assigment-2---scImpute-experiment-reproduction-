# ASSIGNMENT 2 - DATA WRANGLING
# TASK 2 - REVISED, SIMPLIFIED scImpute REPRODUCTION
# Paper: Li & Li (2018), Nature Communications 9:997 (Figure 4a).
# Original algorithm: https://github.com/Vivianstats/scImpute
#
# PURPOSE:
# Re-test scImpute without modifying the authors' algorithm. The first run
# uses UNLABELED cells (the paper's main setting); a second, clearly labeled
# DIAGNOSTIC run provides the true cell types to isolate clustering effects.
#
# IMPORTANT DIFFERENCE FROM THE PAPER:
# The paper multiplied LOG10 expression means by 2..10. Here we apply 2..10
# fold-changes to COUNTS (equivalent to adding log10(2..10) to log10 means).
# This produces more plausible counts without the artificial log10=5 cap
# used in our previous script. We also choose a milder dropout proportion.
# Consequently this is a SIMPLIFIED, ADAPTED experiment, NOT an exact
# numerical reproduction of the published WCSS = 94 / 2646 / 152.
#
# HOW TO RUN:
# 1) Save and extract this .R file in your assignment folder.
# 2) In RStudio, set the working directory to that folder.
# 3) Run the entire file using Source. First use QUICK_TEST <- TRUE.
# 4) If the pilot succeeds, change QUICK_TEST <- FALSE and Source again.
# 5) Upload RESULTS_SUMMARY.csv, CLUSTER_TABLE.csv and FIGURE_PCA.png.
#
# The default uses NEW simulated data. To diagnose ONLY the earlier dataset,
# change DATA_SOURCE to "previous"; the script will reuse the existing count
# matrices from results_scImpute/full/ (no new simulation).


# 1. CHECK / INSTALL THE ORIGINAL scImpute PACKAGE -----------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
user_library <- file.path(path.expand("~"), "R", "assignment2_library")
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_library, .libPaths()))

if (!requireNamespace("scImpute", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = user_library)
  }
  remotes::install_github("Vivianstats/scImpute", lib = user_library,
                          dependencies = NA, upgrade = "never",
                          build_vignettes = FALSE)
}
if (!requireNamespace("scImpute", quietly = TRUE)) {
  stop("The official scImpute package is unavailable. Check installation.")
}
cat("scImpute version:", as.character(packageVersion("scImpute")), "\n")


# 2. DEFINE PARAMETERS AND OUTPUT FOLDERS --------------------------------------
DATA_SOURCE <- "new"             # "new" or "previous" (old input matrices).
QUICK_TEST <- FALSE              # TRUE: 4,000 genes/60 cells; FALSE: 20,000/150.
RUN_LABELED_DIAGNOSTIC <- TRUE  # True labels ONLY in the diagnostic run.
RANDOM_SEED <- 42L
KCLUSTER <- 3L                  # Used only by the unlabeled scImpute run.
DROP_THRESHOLD <- 0.5           # Authors' recommended default.
NCORES <- 1L                    # Safe and reproducible on Windows.
DROP_COEFFICIENT <- 0.1         # From paper: exp(-0.1 * gene_mean^2).
DROPOUT_MULTIPLIER <- 0.70     # OUR choice: milder dropout rate than paper.
PCA_PSEUDOCOUNT <- 1.01

if (!DATA_SOURCE %in% c("new", "previous")) stop("Invalid DATA_SOURCE.")
if (DATA_SOURCE == "previous" && QUICK_TEST) {
  stop("Set QUICK_TEST <- FALSE when DATA_SOURCE = 'previous'.")
}
stopifnot(DROP_THRESHOLD >= 0, DROP_THRESHOLD <= 1,
          DROPOUT_MULTIPLIER > 0, DROPOUT_MULTIPLIER <= 1,
          NCORES >= 1L)

N_TYPES <- 3L
N_GENES <- if (QUICK_TEST) 4000L else 20000L
CELLS_PER_TYPE <- if (QUICK_TEST) 20L else 50L
DE_PER_TYPE <- if (QUICK_TEST) 54L else 270L
N_CELLS <- N_TYPES * CELLS_PER_TYPE

run_type <- if (QUICK_TEST) "pilot" else "full"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_dir <- file.path(getwd(), "results_scImpute_v2",
                      paste0(run_type, "_", DATA_SOURCE, "_", run_stamp))
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
cat("\nWorking directory:", getwd(), "\n")
cat("New experiment results:", base_dir, "\n")


# 3. CREATE (OR REUSE) COMPLETE AND RAW COUNTS ---------------------------------
set.seed(RANDOM_SEED)
cell_types <- rep(paste0("c", seq_len(N_TYPES)), each = CELLS_PER_TYPE)
cell_names <- paste0("Cell_", seq_len(N_CELLS))
gene_names <- paste0("Gene_", seq_len(N_GENES))

if (DATA_SOURCE == "new") {
  # The paper draws gene means and standard deviations on a log10 scale.
  gene_means <- rnorm(N_GENES, mean = 1.8, sd = 0.5)
  gene_sds <- pmax(rnorm(N_GENES, mean = 0.6, sd = 0.1), 0.05)

  # The same number of disjoint DE genes is assigned to each cell type.
  de_groups <- split(sample.int(N_GENES, N_TYPES * DE_PER_TYPE),
                     rep(seq_len(N_TYPES), each = DE_PER_TYPE))

  complete_log10 <- matrix(0, nrow = N_GENES, ncol = N_CELLS,
                           dimnames = list(gene_names, cell_names))
  for (k in seq_len(N_TYPES)) {
    adjusted_mean <- gene_means
    de_genes <- de_groups[[k]]
    # This is an ADAPTATION: multiplicative fold-change on original counts,
    # additive in log space. It avoids unrealistically huge log10 values.
    adjusted_mean[de_genes] <- adjusted_mean[de_genes] +
      log10(sample(2:10, length(de_genes), replace = TRUE))
    inds <- which(cell_types == paste0("c", k))
    complete_log10[, inds] <- matrix(
      rnorm(N_GENES * length(inds),
            mean = rep(adjusted_mean, length(inds)),
            sd = rep(gene_sds, length(inds))),
      nrow = N_GENES, ncol = length(inds))
  }

  # Convert log10(count + 1) to counts, without upper clipping.
  complete_counts <- round(10^complete_log10 - 1)
  complete_counts[complete_counts < 0] <- 0
  dimnames(complete_counts) <- list(gene_names, cell_names)
  # Convert to integer storage only after checking for integer overflow.
  if (max(complete_counts) > .Machine$integer.max) {
    stop("Generated counts exceed R integer range. Inspect the simulation.")
  }
  storage.mode(complete_counts) <- "integer"

  # Each gene/cell pair is independently masked with the gene dropout rate.
  mean_expression <- rowMeans(complete_log10)
  p_dropout <- DROPOUT_MULTIPLIER *
    exp(-DROP_COEFFICIENT * mean_expression^2)
  p_dropout <- pmin(pmax(p_dropout, 0), 1)
  dropout_mask <- matrix(
    runif(N_GENES * N_CELLS) < rep(p_dropout, N_CELLS),
    nrow = N_GENES, ncol = N_CELLS,
    dimnames = dimnames(complete_counts))
  raw_counts <- complete_counts
  raw_counts[dropout_mask] <- 0L

  write.csv(data.frame(gene = gene_names, dropout_probability = p_dropout),
            file.path(base_dir, "SIMULATED_DROPOUT_PROBABILITIES.csv"),
            row.names = FALSE)
  cat("Randomly masked entries:", round(100 * mean(dropout_mask), 2), "%\n")
} else {
  # The previous results remain unchanged. This isolates the effect of
  # using true labels without changing the simulated counts or dropouts.
  previous_dir <- file.path(getwd(), "results_scImpute", "full")
  if (!file.exists(file.path(previous_dir, "complete_counts.csv")) ||
      !file.exists(file.path(previous_dir, "raw_counts.csv"))) {
    stop("Previous count CSV files not found in: ", previous_dir)
  }
  complete_counts <- as.matrix(read.csv(
    file.path(previous_dir, "complete_counts.csv"),
    row.names = 1, check.names = FALSE))
  raw_counts <- as.matrix(read.csv(
    file.path(previous_dir, "raw_counts.csv"),
    row.names = 1, check.names = FALSE))
  storage.mode(complete_counts) <- "numeric"
  storage.mode(raw_counts) <- "numeric"
  gene_names <- rownames(complete_counts)
  cell_names <- colnames(complete_counts)
  if (nrow(complete_counts) != N_GENES || ncol(complete_counts) != N_CELLS) {
    stop("Previous data have unexpected dimensions.")
  }
  # The old experiment used c1=Cell_1..50, c2=51..100, c3=101..150.
  if (!identical(cell_names, paste0("Cell_", seq_len(N_CELLS)))) {
    stop("Previous cell names/order differ. Check the labels before running.")
  }
}

# Verify that all matrices have the correct gene x cell structure.
if (!is.matrix(complete_counts) || !is.matrix(raw_counts) ||
    !identical(dim(complete_counts), c(N_GENES, N_CELLS)) ||
    !identical(dim(raw_counts), dim(complete_counts)) ||
    !identical(dimnames(raw_counts), dimnames(complete_counts))) {
  stop("Invalid gene x cell matrix structure.")
}
stopifnot(all(is.finite(complete_counts)), all(is.finite(raw_counts)),
          all(complete_counts >= 0), all(raw_counts >= 0),
          all(raw_counts <= complete_counts),
          all(colSums(raw_counts) > 0))
removed_positive <- complete_counts > 0 & raw_counts == 0
if (!any(removed_positive)) stop("No positive values were masked.")
cat("Input verified:", N_GENES, "genes x", N_CELLS, "cells\n")

input_file <- file.path(base_dir, "raw_counts.csv")
write.csv(raw_counts, input_file)
write.csv(complete_counts, file.path(base_dir, "complete_counts.csv"))
write.csv(data.frame(cell = cell_names, true_type = cell_types),
          file.path(base_dir, "TRUE_CELL_LABELS.csv"), row.names = FALSE)

# Explicitly check the CSV format expected by the authors' package.
input_check <- read.csv(input_file, row.names = 1, check.names = FALSE)
if (!identical(dim(input_check), c(N_GENES, N_CELLS)) ||
    !identical(rownames(input_check), gene_names) ||
    !identical(colnames(input_check), cell_names)) {
  stop("The saved input CSV does not match the simulation.")
}
rm(input_check)


# 4. DEFINE EVALUATION FUNCTIONS (THE SAME FOR EVERY CONDITION) ---------------
log_transform <- function(x) {
  x <- as.matrix(x)
  if (any(!is.finite(x)) || any(x < 0)) stop("Invalid counts for log10.")
  log10(x + PCA_PSEUDOCOUNT)
}

analyze_counts <- function(counts, labels) {
  # PCA uses CELLS as observations, GENES as variables, independently for
  # each condition. This approximates Figure 4a but is not guaranteed to
  # duplicate all PCA preprocessing decisions from the 2018 paper.
  pca <- prcomp(t(log_transform(counts)), center = TRUE,
                scale. = FALSE, rank. = 2)
  xy <- pca$x[, 1:2, drop = FALSE]
  group_wcss <- vapply(unique(labels), function(label) {
    pts <- xy[labels == label, , drop = FALSE]
    sum(sweep(pts, 2, colMeans(pts), "-")^2)
  }, numeric(1))
  total_ss <- sum(scale(xy, center = TRUE, scale = FALSE)^2)
  wcss <- sum(group_wcss)
  # Normalized within-group fraction reduces sensitivity to the overall
  # PCA scale. Smaller values indicate tighter groups relative to total.
  list(coords = xy, WCSS = wcss, total_ss = total_ss,
       WCSS_fraction = if (total_ss > 0) wcss / total_ss else NA_real_,
       group_wcss = group_wcss)
}

masked_rmse <- function(counts) {
  true_vals <- log_transform(complete_counts)[removed_positive]
  estimated <- log_transform(counts)[removed_positive]
  sqrt(mean((true_vals - estimated)^2))
}

record_metrics <- function(condition, mode, counts, evaluation) {
  data.frame(condition = condition, mode = mode,
             WCSS = evaluation$WCSS,
             WCSS_fraction = evaluation$WCSS_fraction,
             masked_log10_RMSE = masked_rmse(counts),
             zeros_percent = 100 * mean(counts == 0),
             stringsAsFactors = FALSE)
}


# 5. EXECUTE THE AUTHORS' PACKAGE: UNLABELED AND LABELED DIAGNOSTIC ------------
# Crucial distinction: labeled=FALSE is the primary experiment, comparable
# in principle with the paper. labeled=TRUE supplies the true cell types,
# so it is ONLY a diagnostic upper-information scenario, not Figure 4a.
run_scimpute <- function(mode) {
  is_labeled <- identical(mode, "labeled_diagnostic")
  folder <- file.path(base_dir, mode)
  dir.create(folder, showWarnings = FALSE)
  folder <- paste0(normalizePath(folder, winslash = "/", mustWork = TRUE), "/")

  cat("\nRunning scImpute:", mode, "\n")
  set.seed(RANDOM_SEED + if (is_labeled) 102L else 101L)
  if (is_labeled) {
    detected_outliers <- scImpute::scimpute(
      count_path = input_file, infile = "csv", outfile = "csv",
      type = "count", out_dir = folder, labeled = TRUE,
      labels = cell_types, drop_thre = DROP_THRESHOLD, ncores = NCORES)
  } else {
    detected_outliers <- scImpute::scimpute(
      count_path = input_file, infile = "csv", outfile = "csv",
      type = "count", out_dir = folder, labeled = FALSE,
      Kcluster = KCLUSTER, drop_thre = DROP_THRESHOLD, ncores = NCORES)
  }

  result_candidates <- file.path(folder,
    c("scimpute_count.csv", "scImpute_count.csv"))
  matches <- result_candidates[file.exists(result_candidates)]
  if (length(matches) == 0L) {
    stop("scImpute did not save an output CSV in: ", folder)
  }
  counts <- as.matrix(read.csv(matches[1], row.names = 1,
                               check.names = FALSE))
  storage.mode(counts) <- "numeric"
  if (!setequal(rownames(counts), gene_names) ||
      !setequal(colnames(counts), cell_names)) {
    stop("Gene/cell labels in scImpute output do not match input.")
  }
  counts <- counts[gene_names, cell_names, drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Output contains non-finite or negative expression values.")
  }
  write.csv(data.frame(outlier_cell_index = as.integer(detected_outliers)),
            file.path(folder, "OUTLIERS.csv"), row.names = FALSE)

  # The authors' implementation usually saves inferred clusters to clust.rds.
  cluster_file <- file.path(folder, "clust.rds")
  if (file.exists(cluster_file) && !is_labeled) {
    inferred <- readRDS(cluster_file)
    if (length(inferred) == N_CELLS) {
      # Account for named clustering vectors, when supplied by the package.
      if (!is.null(names(inferred)) &&
          setequal(names(inferred), cell_names)) {
        inferred <- inferred[cell_names]
      }
      cluster_table <- as.data.frame.matrix(table(
        True_type = factor(cell_types, levels = unique(cell_types)),
        Inferred_cluster = as.character(inferred), useNA = "ifany"))
      write.csv(cluster_table, file.path(base_dir, "CLUSTER_TABLE.csv"))
      print(table(True_type = cell_types, Inferred_cluster = inferred,
                  useNA = "ifany"))
    } else {
      warning("Unexpected clust.rds length. Inspect it manually.")
    }
  }
  counts
}

# Evaluate complete and raw matrices before any imputation.
evaluation_complete <- analyze_counts(complete_counts, cell_types)
evaluation_raw <- analyze_counts(raw_counts, cell_types)

imputed_unlabeled <- run_scimpute("unlabeled")
evaluation_unlabeled <- analyze_counts(imputed_unlabeled, cell_types)

if (RUN_LABELED_DIAGNOSTIC) {
  imputed_labeled <- run_scimpute("labeled_diagnostic")
  evaluation_labeled <- analyze_counts(imputed_labeled, cell_types)
}


# 6. COMPARE METRICS AND CREATE FIGURES ---------------------------------------
summary_table <- rbind(
  record_metrics("Complete", "reference", complete_counts,
                 evaluation_complete),
  record_metrics("Raw", "reference", raw_counts, evaluation_raw),
  record_metrics("scImpute", "unlabeled", imputed_unlabeled,
                 evaluation_unlabeled))
if (RUN_LABELED_DIAGNOSTIC) {
  summary_table <- rbind(summary_table,
    record_metrics("scImpute", "labeled_diagnostic", imputed_labeled,
                   evaluation_labeled))
}
write.csv(summary_table, file.path(base_dir, "RESULTS_SUMMARY.csv"),
          row.names = FALSE)
cat("\nEXPERIMENTAL RESULTS (not paper values):\n")
print(summary_table, row.names = FALSE)

# Save within-type WCSS so we can identify which population contributes most.
within_table <- data.frame(
  type = unique(cell_types),
  Complete = unname(evaluation_complete$group_wcss),
  Raw = unname(evaluation_raw$group_wcss),
  scImpute_unlabeled = unname(evaluation_unlabeled$group_wcss))
if (RUN_LABELED_DIAGNOSTIC) {
  within_table$scImpute_labeled <- unname(evaluation_labeled$group_wcss)
}
write.csv(within_table, file.path(base_dir, "WITHIN_TYPE_WCSS.csv"),
          row.names = FALSE)

fig_items <- list(
  list(title = "Complete", result = evaluation_complete),
  list(title = "Raw", result = evaluation_raw),
  list(title = "scImpute - unlabeled", result = evaluation_unlabeled))
if (RUN_LABELED_DIAGNOSTIC) {
  fig_items <- append(fig_items, list(
    list(title = "scImpute - labels supplied (diagnostic)",
         result = evaluation_labeled)))
}

all_x <- unlist(lapply(fig_items, function(z) z$result$coords[, 1]))
all_y <- unlist(lapply(fig_items, function(z) z$result$coords[, 2]))
cell_colors <- c(c1 = "#4169A1", c2 = "#26916B", c3 = "#CE6848")

png(file.path(base_dir, "FIGURE_PCA.png"),
    width = if (RUN_LABELED_DIAGNOSTIC) 2000 else 1600,
    height = 580, res = 150)
par(mfrow = c(1, length(fig_items)), mar = c(4.5, 4.3, 3.5, 1))
for (item in fig_items) {
  xy <- item$result$coords
  plot(xy[, 1], xy[, 2], pch = 19, cex = 0.65,
       col = cell_colors[cell_types],
       xlim = range(all_x), ylim = range(all_y),
       main = sprintf("%s\nWCSS = %.1f", item$title, item$result$WCSS),
       xlab = "Principal component 1", ylab = "Principal component 2")
}
par(mfrow = c(1, 1))
dev.off()

# Save coordinates as CSV to allow a more detailed independent inspection.
for (i in seq_along(fig_items)) {
  item <- fig_items[[i]]
  xy <- item$result$coords
  coord_table <- data.frame(cell = cell_names, true_type = cell_types,
                            PC1 = xy[, 1], PC2 = xy[, 2])
  write.csv(coord_table,
            file.path(base_dir, paste0("PCA_COORDINATES_", i, ".csv")),
            row.names = FALSE)
}


# 7. SAVE PARAMETERS, SOFTWARE VERSIONS, AND REPORT NOTES ---------------------
params <- data.frame(
  parameter = c("DATA_SOURCE", "QUICK_TEST", "RUN_LABELED_DIAGNOSTIC",
                "RANDOM_SEED", "N_GENES", "CELLS_PER_TYPE", "DE_PER_TYPE",
                "KCLUSTER", "DROP_THRESHOLD", "NCORES", "DROP_COEFFICIENT",
                "DROPOUT_MULTIPLIER", "PCA_PSEUDOCOUNT", "gene_mean_mu",
                "gene_mean_sd", "gene_sd_mu", "gene_sd_sd", "DE_fold_changes",
                "DE_applied_on", "max_log10_cap", "scImpute_version"),
  value = as.character(c(DATA_SOURCE, QUICK_TEST, RUN_LABELED_DIAGNOSTIC,
                         RANDOM_SEED, N_GENES, CELLS_PER_TYPE, DE_PER_TYPE,
                         KCLUSTER, DROP_THRESHOLD, NCORES, DROP_COEFFICIENT,
                         DROPOUT_MULTIPLIER, PCA_PSEUDOCOUNT, 1.8, 0.5,
                         0.6, 0.1, "2..10", "counts, not log10 means",
                         "none", as.character(packageVersion("scImpute"))))
)
write.csv(params, file.path(base_dir, "PARAMETERS_USED.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(base_dir, "R_SESSION_INFO.txt"))

notes <- c(
  "This is an adapted scImpute simulation, NOT an exact Figure 4a replication.",
  "The original algorithm is unchanged. New input counts avoid arbitrary clipping.",
  "We use count fold-changes (2..10), not multiplication of log10 means.",
  sprintf("Genes: %d; cells: %d; DE genes: %d; random seed: %d.",
          N_GENES, N_CELLS, DE_PER_TYPE * N_TYPES, RANDOM_SEED),
  sprintf("Drop threshold: %.2f; Kcluster: %d; dropout multiplier: %.2f.",
          DROP_THRESHOLD, KCLUSTER, DROPOUT_MULTIPLIER),
  "Unlabeled is the main experiment; labeled is ONLY a diagnostic.",
  "Published Figure 4a: WCSS Complete 94, Raw 2646, scImpute 152.",
  "Do not compare numerical WCSS values as if data and PCA were identical.",
  capture.output(print(summary_table, row.names = FALSE))
)
writeLines(notes, file.path(base_dir, "REPORT_VALUES.txt"))

cat("\nDONE. Main files to send back:\n")
cat("  ", file.path(base_dir, "RESULTS_SUMMARY.csv"), "\n")
cat("  ", file.path(base_dir, "FIGURE_PCA.png"), "\n")
cat("  ", file.path(base_dir, "CLUSTER_TABLE.csv"), "\n")
cat("Also useful: WITHIN_TYPE_WCSS.csv, PARAMETERS_USED.csv, REPORT_VALUES.txt\n")
