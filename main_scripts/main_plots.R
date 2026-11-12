# Generates empirical CDF and bias distribution figures (PDF) from empirical estimation results.
# Set working directory (adjust this path to your own local working directory for replication)
setwd("/Users/gavin/Documents/Antigravity IDE/Thesis/MPI paper/empirical/main_scripts")
library(ggplot2)
library(showtext)
library(reshape2)
library(gridExtra)
library(scales)

USE_RCPP <- TRUE # Set to FALSE to use native R results

font_file <- "../fonts/lmroman10-regular.otf"
if (!file.exists(font_file)) {
  font_file <- "/Library/Fonts/lmroman10-regular.otf"
}
if (file.exists(font_file)) {
  font_add("LM Roman", font_file)
  showtext_auto()
}

csv_file <- if (USE_RCPP) "../output_files/all_households_mpi_rcpp.csv" else "../output_files/all_households_mpi.csv"
df <- read.csv(csv_file)

# 1. MEAN MPI PLOT
plot_data_mean <- data.frame(
  Minimum_MPI = df$Minimum_MPI,
  Maximum_MPI = df$Maximum_MPI,
  Approximated_Mean_MPI = df$Approximated_Mean_MPI,
  Basis_Mean_MPI = df$Basis_Mean_MPI
)

df_melt_mean <- melt(plot_data_mean, variable.name = "MPI_Type", value.name = "MPI_Value")
df_melt_mean <- df_melt_mean[!is.na(df_melt_mean$MPI_Value), ]

p_mean <- ggplot(df_melt_mean, aes(x = MPI_Value, color = MPI_Type)) +
  stat_ecdf(geom = "step", linewidth = 0.5, alpha = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = c("Minimum_MPI" = "blue", "Maximum_MPI" = "red", 
                                "Approximated_Mean_MPI" = "darkgreen", "Basis_Mean_MPI" = "black"),
                     labels = c("Minimum MPI", "Maximum MPI", "Approximate Mean MPI", "Basis Mean MPI")) +
  coord_cartesian(xlim = c(0, 0.4)) +
  labs(x = "Money Pump Index (Mean)", y = "Cumulative Share of Households",
       color = "") +
  theme_minimal(base_family = "LM Roman") +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        panel.grid.minor = element_blank())

dir.create("../tables_figures", showWarnings = FALSE, recursive = TRUE)
pdf("../tables_figures/cdf_plot.pdf", width=4, height=4)
print(p_mean)
dev.off()

# 2. MEDIAN MPI PLOT
plot_data_median <- data.frame(
  Minimum_MPI = df$Minimum_MPI,
  Maximum_MPI = df$Maximum_MPI,
  Approximated_Median_MPI = df$Approximated_Median_MPI,
  Basis_Median_MPI = df$Basis_Median_MPI
)

df_melt_median <- melt(plot_data_median, variable.name = "MPI_Type", value.name = "MPI_Value")
df_melt_median <- df_melt_median[!is.na(df_melt_median$MPI_Value), ]

p_median <- ggplot(df_melt_median, aes(x = MPI_Value, color = MPI_Type)) +
  stat_ecdf(geom = "step", linewidth = 0.5, alpha = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = c("Minimum_MPI" = "blue", "Maximum_MPI" = "red", 
                                "Approximated_Median_MPI" = "darkgreen", "Basis_Median_MPI" = "black"),
                     labels = c("Minimum MPI", "Maximum MPI", "Approximate Median MPI", "Basis Median MPI")) +
  coord_cartesian(xlim = c(0, 0.4)) +
  labs(x = "Money Pump Index (Median)", y = "Cumulative Share of Households",
       color = "") +
  theme_minimal(base_family = "LM Roman") +
  theme(legend.position = c(0.7, 0.2),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        panel.grid.minor = element_blank())

pdf("../tables_figures/cdf_plot_median.pdf", width=4, height=4)
print(p_median)
dev.off()


# 3. EMPIRICAL BIAS HISTOGRAMS
custom_breaks <- function(x) {
  b <- pretty(x)
  if (max(b) < x[2]) b <- c(b, max(b) + diff(b)[1])
  b
}

p1 <- ggplot(df, aes(x = Empirical_Mean_Bias)) + geom_histogram(bins = 100, fill="blue", color="black", alpha=0.7, linewidth=0.2) + 
      theme_minimal(base_family="LM Roman") + theme(plot.title = element_text(hjust = 0.5)) + labs(title="Empirical Mean Bias", x="Bias", y="Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = custom_breaks)
p2 <- ggplot(df, aes(x = Empirical_Median_Bias)) + geom_histogram(bins = 100, fill="red", color="black", alpha=0.7, linewidth=0.2) + 
      theme_minimal(base_family="LM Roman") + theme(plot.title = element_text(hjust = 0.5)) + labs(title="Empirical Median Bias", x="Bias", y="Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = custom_breaks)
p3 <- ggplot(df, aes(x = Mean_Bias_Bound)) + geom_histogram(bins = 100, fill="blue", color="black", alpha=0.4, linewidth=0.2) + 
      theme_minimal(base_family="LM Roman") + theme(plot.title = element_text(hjust = 0.5)) + labs(title="Theoretical Bound for Mean Bias", x="Bound", y="Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = custom_breaks)
p4 <- ggplot(df, aes(x = Median_Bias_Bound)) + geom_histogram(bins = 100, fill="red", color="black", alpha=0.4, linewidth=0.2) + 
      theme_minimal(base_family="LM Roman") + theme(plot.title = element_text(hjust = 0.5)) + labs(title="Theoretical Bound for Median Bias", x="Bound", y="Frequency") + scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = custom_breaks)

pdf("../tables_figures/bias_bounds_hist.pdf", width=8, height=5)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()
