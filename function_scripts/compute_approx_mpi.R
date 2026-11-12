# compute_approx_mpi.R
# Compute approximate MPI over cycles of length 2, 3, and 4 (Echenique et al., 2011)

library(gtools)

# Helper function to generate all possible unique cycles of exactly length 'k'
# from a graph with 'T_val' nodes (observations).
generate_cycles <- function(T_val, k) {
  # If the graph has fewer nodes than the cycle length we want, return empty
  if (T_val < k) return(list())
  
  # Get all unique subsets of nodes of size 'k'
  combos <- combinations(T_val, k)
  cycles <- list()
  
  # For each subset of nodes, generate all valid cyclic permutations
  for(i in 1:nrow(combos)) {
    combo <- combos[i, ]
    fixed <- combo[1] # Fix the first node to avoid counting rotated versions of the same cycle
    rem <- combo[-1]  # Permute the remaining nodes
    perms <- permutations(length(rem), length(rem), rem)
    
    # Store each valid cycle sequence
    for(j in 1:nrow(perms)) {
      cycles[[length(cycles) + 1]] <- c(fixed, perms[j, ])
    }
  }
  return(cycles)
}

library(Rcpp)
    sourceCpp("../function_scripts/rcpp_approx.cpp")

# Evaluates the MPI for all cycles up to length 4 given a cost matrix Z (where Z[i,j] = p_i * q_j).
eval_cycles_for_Z <- function(Z, Q, use_rcpp = FALSE) {
  if (use_rcpp) {
    return(eval_cycles_for_Z_cpp(Z, Q))
  }
  mpis <- numeric(0)
  t_val <- nrow(Z)
  
  # Echenique et al. (2011) approximate the MPI by only looking at cycles of length 2, 3, and 4
  all_cycles <- c(generate_cycles(t_val, 2), 
                  generate_cycles(t_val, 3), 
                  generate_cycles(t_val, 4))
  
  for(c in all_cycles) {
    k <- length(c)
    violation <- TRUE
    cost_sum <- 0
    budget_sum <- 0
    
    # Traverse the edges of the cycle
    for(idx in 1:k) {
      u <- c[idx]
      v <- c[ifelse(idx == k, 1, idx + 1)] # Wrap around to the first node to close the loop
      
      # Check if this edge is a valid revealed preference edge (Z[u,u] >= Z[u,v]). 
      # Also filter zero-cost padding edges (degenerate cycles containing identical bundles)
      # that artificially inflate length.
      if(Z[u, u] < Z[u, v] || all(Q[u,] == Q[v,])) {
        violation <- FALSE
        break
      }
      
      # Accumulate the money pump cost and the total budget across the cycle
      cost_sum <- cost_sum + (Z[u, u] - Z[u, v])
      budget_sum <- budget_sum + Z[u, u]
    }
    
    # If the cycle is a valid GARP violation, calculate its MPI and save it
    if(violation && cost_sum > 1e-7) {
      mpis <- c(mpis, cost_sum / budget_sum)
    }
  }
  
  return(mpis)
}
