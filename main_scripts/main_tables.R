# Generates LaTeX table_1.tex and table_2.tex from empirical estimation results.
# Set working directory (adjust this path to your own local working directory for replication)
setwd("/Users/gavin/Documents/Antigravity IDE/Thesis/MPI paper/empirical/main_scripts")

USE_RCPP <- TRUE # Set to FALSE to use native R results
csv_file <- if (USE_RCPP) "../output_files/all_households_mpi_rcpp.csv" else "../output_files/all_households_mpi.csv"

results_df <- read.csv(csv_file)
cat(sprintf("Loaded %d households from %s\n", nrow(results_df), csv_file))

# Compute aggregate summary statistics across all violating households
res_mat <- as.matrix(results_df[, -1])
stats <- matrix(0, nrow = 8, ncol = 12)
colnames(stats) <- c("Minimum MPI", "Maximum MPI", "Range", 
                     "Approximated Mean MPI", "Approximated Median MPI", 
                     "Basis Mean MPI", "Basis Median MPI", 
                     "Mean Bias Bound", "Median Bias Bound", 
                     "Convergence Iterations", "Empirical Mean Bias", "Empirical Median Bias")
rownames(stats) <- c("Average", "Standard deviation", "Number of zeros", 
                     "Minimum", "First quartile", "Median", "Third quartile", "Maximum")

for(c in 1:12) {
  col_data <- res_mat[, c]
  col_data <- na.omit(col_data)
  
  if(length(col_data) == 0) {
    stats[, c] <- NA
  } else {
    stats[1, c] <- mean(col_data)
    stats[2, c] <- sd(col_data)
    stats[3, c] <- sum(abs(col_data) < 1e-7)
    stats[4, c] <- min(col_data)
    stats[5, c] <- quantile(col_data, 0.25)
    stats[6, c] <- median(col_data)
    stats[7, c] <- quantile(col_data, 0.75)
    stats[8, c] <- max(col_data)
  }
}

# Formatting helpers
fmt <- function(x) {
  x_clean <- ifelse(abs(x) < 0.5e-4, 0, x)
  sprintf("%.4f", x_clean)
}

fmt_row <- function(row_idx, use_dash_except_range=FALSE) {
  vals <- stats[row_idx, ]
  if (use_dash_except_range) {
    return(paste("-", "-", sprintf("%d", as.integer(round(vals[3]))), "-", "-", "-", "-", sep=" & "))
  } else {
    return(paste(fmt(vals[1]), fmt(vals[2]), fmt(vals[3]), fmt(vals[4]), fmt(vals[5]), fmt(vals[6]), fmt(vals[7]), sep=" & "))
  }
}

tab1 <- c(
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}l c c c c c c c}",
  "    \\toprule",
  "    & \\multicolumn{3}{c}{\\textbf{Min/Max MPI}} & \\multicolumn{2}{c}{\\textbf{Approximate MPI}} & \\multicolumn{2}{c}{\\textbf{Basis MPI}} \\\\",
  "    \\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
  "    \\textbf{Statistic} & Min & Max & Range & Average & Median & Average & Median \\\\",
  "    \\midrule",
  paste("    Average &", fmt_row(1), "\\\\"),
  paste("    Standard deviation &", fmt_row(2), "\\\\"),
  paste("    Number of zeros &", fmt_row(3, TRUE), "\\\\"),
  paste("    Minimum &", fmt_row(4), "\\\\"),
  paste("    First quartile &", fmt_row(5), "\\\\"),
  paste("    Median &", fmt_row(6), "\\\\"),
  paste("    Third quartile &", fmt_row(7), "\\\\"),
  paste("    Maximum &", fmt_row(8), "\\\\"),
  "    \\bottomrule",
  "\\end{tabular*}"
)

dir.create("../tables_figures", showWarnings = FALSE, recursive = TRUE)
out_tab1 <- "../tables_figures/table_1.tex"
writeLines(tab1, out_tab1)

fmt_row2 <- function(row_idx, is_count=FALSE) {
  vals <- stats[row_idx, ]
  if (is_count) {
    return(paste(sprintf("%d", as.integer(round(vals[11]))), sprintf("%d", as.integer(round(vals[12]))), sprintf("%d", as.integer(round(vals[8]))), sprintf("%d", as.integer(round(vals[9]))), sep=" & "))
  } else {
    return(paste(fmt(vals[11]), fmt(vals[12]), fmt(vals[8]), fmt(vals[9]), sep=" & "))
  }
}

tab2 <- c(
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}l c c c c}",
  "    \\toprule",
  "    & \\multicolumn{2}{c}{\\textbf{Empirical Bias}} & \\multicolumn{2}{c}{\\textbf{Theoretical Bounds}} \\\\",
  "    \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "    \\textbf{Statistic} & Average & Median & Average & Median \\\\",
  "    \\midrule",
  paste("    Average &", fmt_row2(1), "\\\\"),
  paste("    Standard deviation &", fmt_row2(2), "\\\\"),
  paste("    Number of zeros &", fmt_row2(3, TRUE), "\\\\"),
  paste("    Minimum &", fmt_row2(4), "\\\\"),
  paste("    First quartile &", fmt_row2(5), "\\\\"),
  paste("    Median &", fmt_row2(6), "\\\\"),
  paste("    Third quartile &", fmt_row2(7), "\\\\"),
  paste("    Maximum &", fmt_row2(8), "\\\\"),
  "    \\bottomrule",
  "\\end{tabular*}"
)

out_tab2 <- "../tables_figures/table_2.tex"
writeLines(tab2, out_tab2)

cat(sprintf("Successfully generated %s and %s\n", out_tab1, out_tab2))

