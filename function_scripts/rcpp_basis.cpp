#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

double get_weight(int n, int k) {
    if (k == 2) return 1.0;
    double num = 1.0;
    for (int i = n - k + 2; i <= n; ++i) {
        num *= i;
    }
    return num / k;
}

// [[Rcpp::export]]
NumericVector eval_basis_cycles_for_Z_cpp(NumericMatrix Z, NumericMatrix Q, int max_iter = 5000, int burn_in = 1000) {
    int t_val = Z.nrow();
    if (t_val < 2) return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL);
    
    // Precompute pairwise surpluses, period expenditures, and edge validity
    std::vector<std::vector<double>> surplus(t_val, std::vector<double>(t_val, 0.0));
    std::vector<std::vector<bool>> valid_edge(t_val, std::vector<bool>(t_val, false));
    std::vector<double> exp_vec(t_val, 0.0);
    
    for (int u = 0; u < t_val; ++u) {
        exp_vec[u] = Z(u, u);
        for (int v = 0; v < t_val; ++v) {
            if (u == v) continue;
            surplus[u][v] = Z(u, u) - Z(u, v);
            if (Z(u, u) >= Z(u, v)) {
                bool same_q = true;
                for (int col = 0; col < Q.ncol(); ++col) {
                    if (Q(u, col) != Q(v, col)) {
                        same_q = false;
                        break;
                    }
                }
                valid_edge[u][v] = !same_q;
            }
        }
    }
    
    // Precompute inverse frequency weights by cycle length k
    std::vector<double> weight_by_k(t_val + 1, 0.0);
    for (int k = 2; k <= t_val; ++k) {
        weight_by_k[k] = get_weight(t_val, k);
    }
    
    double cumulative_U = 0;
    double cumulative_W = 0;
    int n_valid = 0;
    int stable_count = 0;
    
    std::vector<double> U_history(max_iter, 0.0);
    std::vector<double> W_history(max_iter, 0.0);
    std::vector<double> rolling_ratios(max_iter, 0.0);
    
    int allocated_size = 50000;
    std::vector<double> all_cycle_mpis;
    all_cycle_mpis.reserve(allocated_size);
    std::vector<double> all_cycle_weights;
    all_cycle_weights.reserve(allocated_size);
    std::vector<int> all_cycle_path_idx;
    all_cycle_path_idx.reserve(allocated_size);
    
    NumericVector src(t_val);
    for (int i = 0; i < t_val; ++i) src[i] = i;
    
    std::vector<double> P_cost(t_val, 0.0);
    std::vector<double> P_exp(t_val + 1, 0.0);
    std::vector<int> P_bad(t_val, 0);
    
    int attempt = 1;
    while (attempt <= max_iter) {
        std::vector<int> perm(t_val);
        NumericVector s = Rcpp::sample(src, t_val, false);
        for (int i = 0; i < t_val; ++i) perm[i] = s[i];
        
        // Precompute cumulative prefix arrays along the sampled Hamiltonian path:
        // P_cost: cumulative pairwise surplus along forward path edges
        // P_exp: cumulative period expenditure along path vertices
        // P_bad: cumulative count of invalid (unaffordable or identical-bundle) forward path edges
        P_cost[0] = 0.0;
        P_bad[0] = 0;
        P_exp[0] = 0.0;
        for (int t = 0; t < t_val; ++t) {
            P_exp[t + 1] = P_exp[t] + exp_vec[perm[t]];
            if (t < t_val - 1) {
                int u = perm[t];
                int v = perm[t + 1];
                P_cost[t + 1] = P_cost[t] + surplus[u][v];
                P_bad[t + 1] = P_bad[t] + (valid_edge[u][v] ? 0 : 1);
            }
        }
        
        double U_p = 0;
        double W_p = 0;
        bool has_valid = false;
        
        // 1. Contiguous subsegment cycles evaluated in O(1) time via prefix sums
        for (int i = 0; i < t_val - 1; ++i) {
            for (int j = i + 1; j < t_val; ++j) {
                int u_last = perm[j];
                int v_first = perm[i];
                if ((P_bad[j] - P_bad[i] == 0) && valid_edge[u_last][v_first]) {
                    double cost_sum = (P_cost[j] - P_cost[i]) + surplus[u_last][v_first];
                    if (cost_sum > 1e-7) {
                        double budget_sum = P_exp[j + 1] - P_exp[i];
                        double mpi_val = cost_sum / budget_sum;
                        int k = j - i + 1;
                        double w = weight_by_k[k];
                        
                        U_p += mpi_val * w;
                        W_p += w;
                        has_valid = true;
                        
                        all_cycle_mpis.push_back(mpi_val);
                        all_cycle_weights.push_back(w);
                        all_cycle_path_idx.push_back(n_valid);
                    }
                }
            }
        }
        
        // 2. Forward skip-edge 2-cycles evaluated in O(1) time
        if (t_val >= 3) {
            double w2 = weight_by_k[2];
            for (int i = 0; i < t_val - 2; ++i) {
                int u = perm[i];
                for (int j = i + 2; j < t_val; ++j) {
                    int v = perm[j];
                    if (valid_edge[u][v] && valid_edge[v][u]) {
                        double cost_sum = surplus[u][v] + surplus[v][u];
                        if (cost_sum > 1e-7) {
                            double budget_sum = exp_vec[u] + exp_vec[v];
                            double mpi_val = cost_sum / budget_sum;
                            
                            U_p += mpi_val * w2;
                            W_p += w2;
                            has_valid = true;
                            
                            all_cycle_mpis.push_back(mpi_val);
                            all_cycle_weights.push_back(w2);
                            all_cycle_path_idx.push_back(n_valid);
                        }
                    }
                }
            }
        }
        
        if (has_valid) {
            n_valid++;
            cumulative_U += U_p;
            cumulative_W += W_p;
            
            U_history[n_valid - 1] = U_p;
            W_history[n_valid - 1] = W_p;
            
            double current_ratio = cumulative_U / cumulative_W;
            rolling_ratios[n_valid - 1] = current_ratio;
            
            if (n_valid > burn_in) {
                double last_ratio = rolling_ratios[n_valid - 2];
                if (std::abs(current_ratio - last_ratio) < 1e-3) {
                    stable_count++;
                } else {
                    stable_count = 0;
                }
                
                if (stable_count >= 200) {
                    break;
                }
            }
        }
        attempt++;
    }
    
    if (n_valid == 0 || cumulative_W == 0.0) return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL);
    
    double mu_U = 0, mu_W = 0;
    for (int i = 0; i < n_valid; ++i) { mu_U += U_history[i]; mu_W += W_history[i]; }
    mu_U /= n_valid; mu_W /= n_valid;
    
    double theta = (mu_W > 0.0) ? (mu_U / mu_W) : 0.0;
    double sig_R = 0.0;
    double sig_W = 0.0;
    double w_min = W_history[0];
    double w_max = W_history[0];
    
    for (int i = 0; i < n_valid; ++i) {
        double R_i = U_history[i] - theta * W_history[i];
        sig_R += R_i * R_i;
        sig_W += (W_history[i] - mu_W) * (W_history[i] - mu_W);
        if (W_history[i] < w_min) w_min = W_history[i];
        if (W_history[i] > w_max) w_max = W_history[i];
    }
    if (n_valid > 1) {
        sig_R = std::sqrt(sig_R / (n_valid - 1));
        sig_W = std::sqrt(sig_W / (n_valid - 1));
    } else {
        sig_R = 0.0;
        sig_W = 0.0;
    }
    
    // Mean bias bound
    double max_bias_bound = (w_min > 0.0 && mu_W > 0.0 && n_valid > 0) ? 
        ((sig_R * sig_W) / (mu_W * w_min * n_valid)) : 0.0;
    if (std::isnan(max_bias_bound) || std::isinf(max_bias_bound)) max_bias_bound = 0.0;
    max_bias_bound = std::min(1.0, max_bias_bound);
    
    double final_ratio = cumulative_U / cumulative_W;
    
    // Sort extracted cycles by MPI
    int total_cycles = all_cycle_mpis.size();
    std::vector<int> sort_idx(total_cycles);
    std::iota(sort_idx.begin(), sort_idx.end(), 0);
    std::sort(sort_idx.begin(), sort_idx.end(), [&](int a, int b) {
        return all_cycle_mpis[a] < all_cycle_mpis[b];
    });
    
    double total_weight = 0;
    for (double w : all_cycle_weights) total_weight += w;
    double target_weight = total_weight / 2.0;
    
    // Midpoint weighted median
    double m_lower = 0.0;
    double m_upper = 0.0;
    double cum_weight = 0.0;
    bool found_lower = false;
    
    for (int i = 0; i < total_cycles; ++i) {
        cum_weight += all_cycle_weights[sort_idx[i]];
        if (cum_weight >= target_weight && !found_lower) {
            m_lower = all_cycle_mpis[sort_idx[i]];
            found_lower = true;
            if (cum_weight > target_weight || i == total_cycles - 1) {
                m_upper = m_lower;
                break;
            } else {
                m_upper = all_cycle_mpis[sort_idx[i + 1]];
                break;
            }
        }
    }
    double weighted_median = 0.5 * (m_lower + m_upper);
    
    // Sample delta gap for median (alpha = 0.50)
    double delta_val = 0.0;
    if (total_weight > 0.0) {
        double weight_le = 0.0;
        double weight_lt = 0.0;
        for (int i = 0; i < total_cycles; ++i) {
            double mpi_i = all_cycle_mpis[i];
            double w_i = all_cycle_weights[i];
            if (mpi_i <= weighted_median) weight_le += w_i;
            if (mpi_i < weighted_median - 1e-12) weight_lt += w_i;
        }
        double cdf_at_median = weight_le / total_weight;
        double cdf_prev = weight_lt / total_weight;
        double upper_gap = cdf_at_median - 0.50;
        double lower_gap = 0.50 - cdf_prev;
        if (upper_gap > 0.0 && lower_gap > 0.0) {
            delta_val = std::min(upper_gap, lower_gap);
        }
    }
    
    // Median bias bound
    double median_bias_bound = 1.0;
    if (delta_val > 0.0 && w_max > 0.0) {
        median_bias_bound = 4.0 * std::exp(- (2.0 * n_valid * delta_val * delta_val * w_min * w_min) / (w_max * w_max));
        median_bias_bound = std::min(1.0, median_bias_bound);
    }
    
    // Bootstrap bias estimation
    int B = 1000;
    double sum_boot_means = 0.0;
    double sum_boot_medians = 0.0;
    
    NumericVector src_paths(n_valid);
    for (int i = 0; i < n_valid; ++i) src_paths[i] = i;
    std::vector<int> path_counts(n_valid, 0);
    
    for (int b = 0; b < B; ++b) {
        // Resample H-path bases with replacement
        NumericVector b_idx_paths = Rcpp::sample(src_paths, n_valid, true);
        std::fill(path_counts.begin(), path_counts.end(), 0);
        double b_sum_U = 0.0;
        double b_sum_W = 0.0;
        for (int i = 0; i < n_valid; ++i) {
            int idx = b_idx_paths[i];
            path_counts[idx]++;
            b_sum_U += U_history[idx];
            b_sum_W += W_history[idx];
        }
        sum_boot_means += (b_sum_U / b_sum_W);
        
        // Midpoint weighted median over cycles in the resampled H-path bases
        double b_med = weighted_median;
        if (b_sum_W > 0.0) {
            double b_target = b_sum_W / 2.0;
            double b_cum = 0.0;
            double b_lower = weighted_median;
            double b_upper = weighted_median;
            for (int i = 0; i < total_cycles; ++i) {
                int c = sort_idx[i];
                int mult = path_counts[all_cycle_path_idx[c]];
                if (mult == 0) continue;
                b_cum += all_cycle_weights[c] * mult;
                if (b_cum >= b_target) {
                    b_lower = all_cycle_mpis[c];
                    if (b_cum > b_target) {
                        b_upper = b_lower;
                    } else {
                        b_upper = b_lower;
                        for (int j = i + 1; j < total_cycles; ++j) {
                            int c_next = sort_idx[j];
                            if (path_counts[all_cycle_path_idx[c_next]] > 0) {
                                b_upper = all_cycle_mpis[c_next];
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            b_med = 0.5 * (b_lower + b_upper);
        }
        sum_boot_medians += b_med;
    }
    
    double empirical_mean_bias = (sum_boot_means / B) - final_ratio;
    double empirical_median_bias = (sum_boot_medians / B) - weighted_median;
    
    int iterations_completed = (attempt <= max_iter) ? attempt : max_iter;
    return NumericVector::create(final_ratio, weighted_median, max_bias_bound, median_bias_bound, iterations_completed, empirical_mean_bias, empirical_median_bias);
}
