#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Helper to evaluate cycles of a specific length k
std::vector<double> generate_and_eval_cycles(int t_val, int k, const NumericMatrix& Z, const NumericMatrix& Q) {
    std::vector<double> mpis;
    if (t_val < k) return mpis;
    
    // We want all combinations of k nodes from 1 to t_val
    // A standard way is to use a boolean mask of length t_val with k trues,
    // then permute the mask to generate combinations.
    std::string bitmask(k, 1);
    bitmask.resize(t_val, 0);
    
    do {
        std::vector<int> combo;
        for (int i = 0; i < t_val; ++i) {
            if (bitmask[i]) combo.push_back(i);
        }
        
        // For each combination, we fix the first element and permute the rest
        // to generate all unique cyclic permutations.
        int fixed = combo[0];
        std::vector<int> rem(combo.begin() + 1, combo.end());
        
        // Sort to ensure we generate all permutations
        std::sort(rem.begin(), rem.end());
        
        do {
            std::vector<int> cycle;
            cycle.push_back(fixed);
            for(int val : rem) cycle.push_back(val);
            
            bool violation = true;
            double cost_sum = 0;
            double budget_sum = 0;
            
            for(int idx = 0; idx < k; ++idx) {
                int u = cycle[idx];
                int v = cycle[(idx == k - 1) ? 0 : idx + 1];
                
                bool identical = true;
                for (int col = 0; col < Q.ncol(); ++col) {
                    if (Q(u, col) != Q(v, col)) {
                        identical = false;
                        break;
                    }
                }
                
                if (Z(u, u) < Z(u, v) || identical) {
                    violation = false;
                    break;
                }
                
                cost_sum += (Z(u, u) - Z(u, v));
                budget_sum += Z(u, u);
            }
            
            if (violation && cost_sum > 1e-7) {
                mpis.push_back(cost_sum / budget_sum);
            }
            
        } while (std::next_permutation(rem.begin(), rem.end()));
        
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    
    return mpis;
}

// [[Rcpp::export]]
NumericVector eval_cycles_for_Z_cpp(NumericMatrix Z, NumericMatrix Q) {
    int t_val = Z.nrow();
    std::vector<double> all_mpis;
    
    for (int k = 2; k <= 4; ++k) {
        std::vector<double> mpis_k = generate_and_eval_cycles(t_val, k, Z, Q);
        all_mpis.insert(all_mpis.end(), mpis_k.begin(), mpis_k.end());
    }
    
    return wrap(all_mpis);
}
