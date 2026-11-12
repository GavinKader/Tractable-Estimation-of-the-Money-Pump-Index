# Main script to replicate the empirical analysis in "Tractable Estimation of the Money Pump Index".
# Set working directory (adjust this path to your own local working directory for replication)
setwd("/Users/gavin/Documents/Antigravity IDE/Thesis/MPI paper/empirical/main_scripts")
library(parallel)
source("../function_scripts/compute_approx_mpi.R")
source("../function_scripts/compute_minmax_mpi.R")
source("../function_scripts/compute_basis_mpi.R")

USE_RCPP <- TRUE # Set to FALSE to run in native R

# Load Data
price_raw <- read.table("../data/price.raw", header = FALSE)
purchase_raw <- read.table("../data/purchase.raw", header = FALSE)

T_val <- 26 # Observations
N_val <- 375 # Goods
NumPanelists <- 494 # Households

# Build the global price matrix P (Observations x Goods)
P_matrix <- matrix(0, nrow=T_val, ncol=N_val)
for(i in 1:nrow(price_raw)) {
  item <- price_raw[i, 1]
  week <- price_raw[i, 2]
  price <- price_raw[i, 3]
  P_matrix[week, item] <- price
}

# Build a list of quantity matrices Q for each household
Q_list <- vector("list", NumPanelists)
for(i in 1:NumPanelists) {
  Q_list[[i]] <- matrix(0, nrow=T_val, ncol=N_val)
}

for(i in 1:nrow(purchase_raw)) {
  item <- purchase_raw[i, 1]
  week <- purchase_raw[i, 2]
  qty <- purchase_raw[i, 3]
  panelist <- purchase_raw[i, 4]
  Q_list[[panelist]][week, item] <- qty
}

# Pre-compute valid Z matrices (Cost Matrices) for all panelists.
# Z[i,j] represents the cost of buying the bundle from week 'j' at week 'i's prices and filtering out observations with zero total consumption.
Z_list <- list()
Q_valid_list <- list()

for(i in 1:NumPanelists) {
  Q <- Q_list[[i]]
  keep <- rowSums(Q) > 0
  p_new <- P_matrix[keep, , drop = FALSE]
  q_new <- Q[keep, , drop = FALSE]
  
  if (nrow(p_new) > 1) {
    Z_list[[i]] <- p_new %*% t(q_new)
    Q_valid_list[[i]] <- q_new
  } else {
    Z_list[[i]] <- NULL
    Q_valid_list[[i]] <- NULL
  }
}

# Enforce parallel reproducibility via L'Ecuyer-CMRG
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
num_cores <- detectCores() - 1 

# PHASE 1: MIN AND MAX MPI
cat("--- Phase 1: Min and Max MPI ---\n")
start_minmax <- Sys.time()

process_minmax <- function(i) {
  Z <- Z_list[[i]]
  Q_valid <- Q_valid_list[[i]]
  if (is.null(Z)) return(NULL)
  
  mm <- compute_min_max_mpi(Z, Q_valid, use_rcpp = USE_RCPP)
  if(mm[1] == 0 && mm[2] == 0) return(NULL) # Perfectly rational
  
  return(list(panelist = i, min_mpi = mm[1], max_mpi = mm[2], range = mm[2] - mm[1]))
}

res_minmax <- mclapply(1:NumPanelists, process_minmax, mc.cores = num_cores)
violating_minmax <- res_minmax[!sapply(res_minmax, is.null)]
violating_ids <- sapply(violating_minmax, function(x) x$panelist)

end_minmax <- Sys.time()
cat(sprintf("Time taken for Min/Max MPI: %.2f seconds\n\n", 
            as.numeric(difftime(end_minmax, start_minmax, units="secs"))))

# Phase 2: Approximated MPI (cycles up to length 4)
cat("--- Phase 2: Approximated Mean and Median MPI (All Cycles up to length 4) ---\n")
start_approx <- Sys.time()

process_approx <- function(i) {
  Z <- Z_list[[i]]
  Q_valid <- Q_valid_list[[i]]
  mpis <- eval_cycles_for_Z(Z, Q_valid, use_rcpp = USE_RCPP)
  
  if(length(mpis) == 0) {
    return(c(NA, NA))
  } else {
    return(c(mean(mpis), median(mpis)))
  }
}

# Only process violating panelists
res_approx <- mclapply(violating_ids, process_approx, mc.cores = num_cores)

end_approx <- Sys.time()
cat(sprintf("Time taken for Approximated MPI: %.2f seconds\n\n", 
            as.numeric(difftime(end_approx, start_approx, units="secs"))))

# Phase 3: Basis MPI (sample Hamiltonian paths)
cat("--- Phase 3: Basis Mean and Median MPI (Fundamental Cycle Basis) ---\n")
start_basis <- Sys.time()

process_basis <- function(i) {
  Z <- Z_list[[i]]
  Q_valid <- Q_valid_list[[i]]
  # Set iteration and burn-in parameters
  res <- eval_basis_cycles_for_Z(Z, Q_valid, max_iter = 5000, burn_in = 1000, use_rcpp = USE_RCPP)
  return(res)
}

# Only process violating panelists
res_basis <- mclapply(violating_ids, process_basis, mc.cores = num_cores)

end_basis <- Sys.time()
cat(sprintf("Time taken for Basis MPI: %.2f seconds\n\n", 
            as.numeric(difftime(end_basis, start_basis, units="secs"))))

# COMPILE RESULTS AND GENERATE TABLE
res_mat <- matrix(0, nrow = length(violating_ids), ncol = 12)
for (idx in 1:length(violating_ids)) {
  res_mat[idx, 1] <- violating_minmax[[idx]]$min_mpi
  res_mat[idx, 2] <- violating_minmax[[idx]]$max_mpi
  res_mat[idx, 3] <- violating_minmax[[idx]]$range
  res_mat[idx, 4] <- res_approx[[idx]][1]
  res_mat[idx, 5] <- res_approx[[idx]][2]
  res_mat[idx, 6] <- res_basis[[idx]][1]
  res_mat[idx, 7] <- res_basis[[idx]][2]
  res_mat[idx, 8] <- res_basis[[idx]][3] # Mean Bias Bound
  res_mat[idx, 9] <- res_basis[[idx]][4] # Median Bias Bound
  res_mat[idx, 10] <- res_basis[[idx]][5] # Convergence Iterations
  res_mat[idx, 11] <- res_basis[[idx]][6] # Empirical Mean Bias (Bootstrap)
  res_mat[idx, 12] <- res_basis[[idx]][7] # Empirical Median Bias (Bootstrap)
}

cat(sprintf("Done! Processed %d violating panelists.\n\n", length(violating_ids)))

# Save the raw matrix of all household MPIs to a CSV for caching and future analysis
results_df <- data.frame(
  Household = violating_ids,
  Minimum_MPI = numeric(length(violating_ids)),
  Maximum_MPI = numeric(length(violating_ids)),
  Range = numeric(length(violating_ids)),
  Approximated_Mean_MPI = numeric(length(violating_ids)),
  Approximated_Median_MPI = numeric(length(violating_ids)),
  Basis_Mean_MPI = numeric(length(violating_ids)),
  Basis_Median_MPI = numeric(length(violating_ids)),
  Mean_Bias_Bound = numeric(length(violating_ids)),
  Median_Bias_Bound = numeric(length(violating_ids)),
  Convergence_Iterations = numeric(length(violating_ids)),
  Empirical_Mean_Bias = numeric(length(violating_ids)),
  Empirical_Median_Bias = numeric(length(violating_ids))
)

for (idx in 1:length(violating_ids)) {
  results_df[idx, "Minimum_MPI"] <- violating_minmax[[idx]]$min_mpi
  results_df[idx, "Maximum_MPI"] <- violating_minmax[[idx]]$max_mpi
  results_df[idx, "Range"] <- violating_minmax[[idx]]$range
  results_df[idx, "Approximated_Mean_MPI"] <- res_approx[[idx]][1]
  results_df[idx, "Approximated_Median_MPI"] <- res_approx[[idx]][2]
  results_df[idx, "Basis_Mean_MPI"] <- res_basis[[idx]][1]
  results_df[idx, "Basis_Median_MPI"] <- res_basis[[idx]][2]
  results_df[idx, "Mean_Bias_Bound"] <- res_basis[[idx]][3]
  results_df[idx, "Median_Bias_Bound"] <- res_basis[[idx]][4]
  results_df[idx, "Convergence_Iterations"] <- res_basis[[idx]][5]
  results_df[idx, "Empirical_Mean_Bias"] <- res_basis[[idx]][6]
  results_df[idx, "Empirical_Median_Bias"] <- res_basis[[idx]][7]
}

out_csv <- if(USE_RCPP) "../output_files/all_households_mpi_rcpp.csv" else "../output_files/all_households_mpi.csv"
write.csv(results_df, out_csv, row.names = FALSE)
cat("Individual household MPI estimates saved to ", out_csv, "\n\n")

# Compute aggregate summary statistics across all violating households
res_mat <- as.matrix(results_df[, -1])
stats <- matrix(0, nrow=8, ncol=12)
colnames(stats) <- c("Minimum MPI", "Maximum MPI", "Range", "Approximated Mean MPI", "Approximated Median MPI", "Basis Mean MPI", "Basis Median MPI", "Mean Bias Bound", "Median Bias Bound", "Convergence Iterations", "Empirical Mean Bias", "Empirical Median Bias")
rownames(stats) <- c("Average", "Standard deviation", "Number of zeros", "Minimum", "First quartile", "Median", "Third quartile", "Maximum")

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

out_table <- if(USE_RCPP) "../output_files/table1_replication_rcpp.csv" else "../output_files/table1_replication.csv"
write.csv(round(stats, 4), out_table)

time_minmax <- as.numeric(difftime(end_minmax, start_minmax, units="secs"))
time_approx <- as.numeric(difftime(end_approx, start_approx, units="secs"))
time_basis <- as.numeric(difftime(end_basis, start_basis, units="secs"))
times_df <- data.frame(Phase = c("MinMax", "Approx", "Basis"), Time_Seconds = c(time_minmax, time_approx, time_basis))
out_times <- if(USE_RCPP) "../output_files/times_rcpp.csv" else "../output_files/times_r.csv"
write.csv(times_df, out_times, row.names = FALSE)

cat("\nNote 1: Run 'Rscript main_tables.R' to generate table_1.tex and table_2.tex\n")
cat("\nNote 2: Both native R and rcpp versions of the functions create identical outputs for the deterministic calculations.\n")
cat("Run 'Rscript main_plots.R' to generate empirical CDF and bias figures.\n\n")

# Phase 4: Validation of exact MPIs
# We verify if the estimator correctly collapses to this exact value without bias for the cases where min and max bound coincide.
cat("--- Phase 4: Exact Theoretical Validation ---\n")
exact_mask <- abs(res_mat[, 1] - res_mat[, 2]) < 1e-7
exact_count <- sum(exact_mask, na.rm = TRUE)

cat(sprintf("Found %d households where Min MPI == Max MPI.\n", exact_count))

if(exact_count > 0) {
  exact_min_mpis <- res_mat[exact_mask, 1]
  exact_basis_means <- res_mat[exact_mask, 6]
  exact_basis_medians <- res_mat[exact_mask, 7]
  
  for(i in 1:exact_count) {
    cat(sprintf("Household %d: Min = %.5f, Basis Mean = %.5f, Basis Median = %.5f\n", 
                i, exact_min_mpis[i], exact_basis_means[i], exact_basis_medians[i]))
  }
  
  max_dev_mean <- max(abs(exact_min_mpis - exact_basis_means), na.rm = TRUE)
  max_dev_med <- max(abs(exact_min_mpis - exact_basis_medians), na.rm = TRUE)
  mse_mean <- mean((exact_min_mpis - exact_basis_means)^2, na.rm = TRUE)
  mse_med <- mean((exact_min_mpis - exact_basis_medians)^2, na.rm = TRUE)
  
  cat(sprintf("Validation Results against Ground-Truth Exact MPI:\n"))
  cat(sprintf("Basis Mean - Maximum Absolute Deviation: %.10f, MSE: %.10f\n", max_dev_mean, mse_mean))
  cat(sprintf("Basis Median - Maximum Absolute Deviation: %.10f, MSE: %.10f\n", max_dev_med, mse_med))
  
  identical_idx <- which(results_df$Range < 1e-7)
  cat(sprintf("Households with identical Min/Max MPI: %d\n", length(identical_idx)))
  cat(sprintf("Households where Basis Mean matches Min/Max exactly: %d\n", 
              sum(abs(results_df$Basis_Mean_MPI[identical_idx] - results_df$Minimum_MPI[identical_idx]) < 1e-5)))
  cat(sprintf("Households where Basis Median matches Min/Max exactly: %d\n", 
              sum(abs(results_df$Basis_Median_MPI[identical_idx] - results_df$Minimum_MPI[identical_idx]) < 1e-5)))
  
  if (max_dev_mean < 1e-3 && max_dev_med < 1e-3) {
    cat("100% match for both Basis Mean and Basis Median!\n")
  } else {
    cat("Deviation detected\n")
  }
}