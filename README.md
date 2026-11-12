# Tractable Estimation of the Money Pump Index

Replication package for the empirical analysis, tables, and figures in *"Tractable Estimation of the Money Pump Index"* ([arXiv:2609.16086](https://arxiv.org/abs/2609.16086)) by Gavin Kader/

## Repository Structure

```text
.
├── data/
│   ├── price.raw                    # Stanford Basket item prices (26 weeks × 375 goods)
│   └── purchase.raw                 # Stanford Basket household purchase quantities (494 households)
├── fonts/
│   └── lmroman10-regular.otf        # Latin Modern Roman font used for plot typography
├── function_scripts/
│   ├── compute_minmax_mpi.R         # R wrapper & native R Min/Max MPI bounds (Smeulders et al., 2013)
│   ├── rcpp_minmax.cpp              # Rcpp C++ backend for Min/Max MPI bounds
│   ├── compute_approx_mpi.R         # R wrapper & native R short-cycle MPI (k ≤ 4; Echenique et al., 2011)
│   ├── rcpp_approx.cpp              # Rcpp C++ backend for short-cycle approximated MPI
│   ├── compute_basis_mpi.R          # R wrapper & native R Basis Mean/Median MPI + bias bounds & bootstrap
│   └── rcpp_basis.cpp               # Rcpp C++ backend for Basis Mean/Median MPI + bias bounds & bootstrap
├── main_scripts/
│   ├── main.R                       # Runs Min/Max, Approximated, and Basis MPI across all households
│   ├── main_tables.R                # Generates LaTeX tables (table_1.tex, table_2.tex)
│   └── main_plots.R                 # Generates empirical CDF and bias distribution plots (PDF)
├── output_files/
│   ├── all_households_mpi_rcpp.csv  # Per-household MPI estimates & bias metrics (Rcpp run, N = 396)
│   ├── all_households_mpi.csv       # Per-household MPI estimates & bias metrics (native R run, N = 396)
│   ├── table1_replication_rcpp.csv  # Summary statistics matrix (Rcpp run)
│   ├── table1_replication.csv       # Summary statistics matrix (native R run)
│   ├── times_rcpp.csv               # Phase execution times (Rcpp run)
│   └── times_r.csv                  # Phase execution times (native R run)
└── tables_figures/
    ├── table_1.tex                  # Generated LaTeX Table 1 (MPI descriptive statistics)
    ├── table_2.tex                  # Generated LaTeX Table 2 (empirical bias & theoretical bounds)
    ├── cdf_plot.pdf                 # Empirical CDF of Mean MPI across GARP-violating households
    ├── cdf_plot_median.pdf          # Empirical CDF of Median MPI across GARP-violating households
    └── bias_bounds_hist.pdf         # Histograms of empirical bootstrap bias & analytical bias bounds
```

## Replication

### 1. Requirements

- **R:** Tested on R version `4.5.3` (`2026-03-11`, `aarch64-apple-darwin20`).
- **C++ compiler toolchain:** Required for `Rcpp::sourceCpp()` (e.g., Xcode Command Line Tools on macOS, `build-essential` on Linux, or Rtools on Windows). Because the scripts in `function_scripts/` compile the `.cpp` files via `sourceCpp()` when sourced, a working C++ compiler is needed.
- **R packages:**
  - `main.R` & `function_scripts/`: `Rcpp`, `gtools`, `parallel` (base R)
  - `main_plots.R`: `ggplot2`, `showtext`, `reshape2`, `gridExtra`, `scales`
  - `main_tables.R`: base R only

Install the CRAN packages in R with:

```r
install.packages(c("Rcpp", "gtools", "ggplot2", "showtext", "reshape2", "gridExtra", "scales"))
```

> **OS Note:** `main.R` uses forked parallelism (`parallel::mclapply` with `mc.cores = detectCores() - 1`) and `RNGkind("L'Ecuyer-CMRG")` for parallel reproducibility on macOS and Linux. On Windows, `mclapply` does not support `mc.cores > 1`; set `num_cores <- 1` in `main_scripts/main.R` before running.

### 2. Running the Replication Scripts

1. **Set the working directory:** Update the `setwd(...)` path on line 3 of `main_scripts/main.R`, `main_scripts/main_tables.R`, and `main_scripts/main_plots.R` to point to your local `main_scripts/` directory (or remove/comment out the `setwd(...)` line and execute from inside `main_scripts/`).
2. **Select the computation backend (`USE_RCPP`):**
   - `main_scripts/main.R`, `main_scripts/main_tables.R`, and `main_scripts/main_plots.R` all default to `USE_RCPP <- TRUE` (`TRUE` runs compiled C++ via `Rcpp` and reads/writes `*_rcpp.csv`; set `USE_RCPP <- FALSE` to run or load the native-R backend instead).
   - Both backends yield identical deterministic values for Min/Max MPI and Approximated MPI (`k ≤ 4`). For Basis MPI, the Rcpp backend draws Hamiltonian paths from the `L'Ecuyer-CMRG` seed (`12345`) set in `main.R`, whereas the native-R `eval_basis_cycles_for_Z()` function sets `set.seed(2026)` internally.
3. **Execute the scripts in order:**
   1. `main_scripts/main.R`
   2. `main_scripts/main_tables.R`
   3. `main_scripts/main_plots.R`

### 3. Outputs

All replication outputs are self-contained within the repository:

- **`main.R`** processes all `494` households across `26` weeks and `375` goods, identifies the `396` GARP-violating households, validates agreement between the Basis estimates and exact MPI for the `75` households where Minimum and Maximum MPI coincide, and writes results to `../output_files/`.
- **`main_tables.R`** writes `table_1.tex` and `table_2.tex` to `../tables_figures/`.
- **`main_plots.R`** writes `cdf_plot.pdf`, `cdf_plot_median.pdf`, and `bias_bounds_hist.pdf` to `../tables_figures/`.

### 4. Logged Reference Runtimes

Reference execution times recorded in `output_files/times_rcpp.csv` and `output_files/times_r.csv`:
## Data Source

The files `data/price.raw` and `data/purchase.raw` contain the Stanford Basket scanner dataset originally analyzed in Echenique, Lee, and Shum (2011) and subsequently in Smeulders et al. (2013), covering weekly grocery expenditures for 494 households over 26 weeks and 375 goods.

## References

- Echenique, F., Lee, S., & Shum, M. (2011). The money pump as a measure of revealed preference violations. *Journal of Political Economy*, 119(6), 1201–1223.

- Smeulders, B., Cherchye, L., De Rock, B., & Spieksma, F. C. R. (2013). The money pump as a measure of revealed preference violations: A comment. *Journal of Political Economy*, 121(6), 1248–1258.

- Kader, G. (2026). *Tractable Estimation of the Money Pump Index*. arXiv preprint [arXiv:2609.16086](https://arxiv.org/abs/2609.16086).
