# compute_minmax_mpi.R
# Compute minimum and maximum MPI via Bellman-Ford binary search (Smeulders et al., 2013)

# Helper function to compute the indegree (number of incoming edges) for each node.
# Used for topological sorting to detect cycles.
indegree <- function(g) {
  incoming <- numeric(length(g))
  for (k in 1:length(g)) {
    if(length(g[[k]]) > 0) {
      for(j in 1:length(g[[k]])) {
        i <- g[[k]][[j]][1]
        incoming[i] <- incoming[i] + 1
      }
    }
  }
  return(incoming)
}

# Determines if the graph has any cycles by attempting a topological sort.
# If we can topologically sort all nodes, the graph is a Directed Acyclic Graph (DAG).
# Returns 1 if a cycle exists, 0 if acyclic.
cycle_detection_topo <- function(Graph) {
  incoming <- indegree(Graph)
  lst <- which(incoming == 0) # Start with nodes that have no incoming edges
  
  count <- 0
  while(length(lst) > 0) {
    node <- lst[length(lst)]
    lst <- lst[-length(lst)]
    count <- count + 1
    
    # Remove this node and decrement indegree of its neighbors
    if(length(Graph[[node]]) > 0) {
      for(j in 1:length(Graph[[node]])) {
        i <- Graph[[node]][[j]][1]
        incoming[i] <- incoming[i] - 1
        if(incoming[i] == 0) {
          lst <- c(lst, i)
        }
      }
    }
  }
  
  # If we couldn't process all nodes, there must be a cycle trapping the remaining nodes
  if(count < length(Graph)) return(1)
  else return(0)
}

# Standard Bellman-Ford algorithm to compute shortest paths and detect negative weight cycles.
# Returns -1 immediately if it detects a negative weight cycle.
bellmanford <- function(g, s) {
  n <- length(g)
  d <- rep(1e6, n) # Initialize distances to infinity
  pred <- rep(0, n)
  d[s] <- 0
  
  # Relax edges (n-1) times
  for(t in 1:(n-1)) {
    for(i in 1:n) {
      if(length(g[[i]]) > 0) {
        for(k in 1:length(g[[i]])) {
          j <- g[[i]][[k]][1]
          w <- g[[i]][[k]][2]
          if(d[i] + w < d[j]) {
            d[j] <- d[i] + w
            pred[j] <- i
          }
        }
      }
    }
  }
  
  # Check for negative weight cycles
  for(i in 1:n) {
    if(length(g[[i]]) > 0) {
      for(k in 1:length(g[[i]])) {
        j <- g[[i]][[k]][1]
        w <- g[[i]][[k]][2]
        if(d[i] + w < d[j]) {
          return(-1) # Negative cycle found
        }
      }
    }
  }
  
  return(list(d=d, pred=pred))
}

# Constructs the zero-residual graph from the shortest path distances.
# Edges with a reduced cost of exactly zero are kept.
zero_residual_graph <- function(Graph, d) {
  g <- vector("list", length(Graph))
  for(i in 1:length(Graph)) {
    g[[i]] <- list()
    if(length(Graph[[i]]) > 0) {
      for(k in 1:length(Graph[[i]])) {
        j <- Graph[[i]][[k]][1]
        w <- Graph[[i]][[k]][2]
        if(abs(w + d[i] - d[j]) < 1e-9) {
          g[[i]][[length(g[[i]]) + 1]] <- c(j, 0)
        }
      }
    }
  }
  return(g)
}

# Adjust edge weights by subtracting u * budget
redefined_graph <- function(Graph, u) {
  g <- Graph
  for(i in 1:length(g)) {
    if(length(g[[i]]) > 0) {
      for(k in 1:length(g[[i]])) {
        g[[i]][[k]][2] <- g[[i]][[k]][2] - u * g[[i]][[k]][3]
      }
    }
  }
  return(g)
}

# Uses binary search to find the precise ratio 'u' (the MPI) 
# that drives the shortest-path graph to transition.
minimum_cost_time <- function(Graph, low, high) {
  if(cycle_detection_topo(Graph) == 0) return(0)
  
  while(abs(high - low) > 1e-10) {
    u <- (low + high) / 2
    g <- redefined_graph(Graph, u)
    res <- bellmanford(g, 1)
    
    if(is.list(res)) {
      d <- res$d
      e <- zero_residual_graph(g, d)
      if(cycle_detection_topo(e) == 1) {
        return(u)
      } else {
        low <- u
      }
    } else {
      high <- u
    }
  }
  return((high + low) / 2)
}

library(Rcpp)
sourceCpp("../function_scripts/rcpp_minmax.cpp")

# Main wrapper to compute both the minimum and maximum MPI for a household.
# It formulates the problem in two different ways to run the binary search for both bounds.
compute_min_max_mpi <- function(Z, Q, use_rcpp = FALSE) {
  if (use_rcpp) {
    return(compute_min_max_mpi_cpp(Z, Q))
  }
  t <- nrow(Z)
  G <- matrix(0, nrow=t, ncol=t)
  
  # Base strict revealed preference graph
  for(i in 1:t) {
    for(j in 1:t) {
      if(Z[i,j] <= Z[i,i]) G[i,j] <- 1
    }
  }
  
  # Graph for Minimum MPI
  g_min <- vector("list", t)
  for(i in 1:t) {
    for(k in 1:t) {
      # Add edge if strict preference holds, but filter out padding edges (identical bundles)
      if(k != i && G[i,k] == 1 && !all(Q[i,] == Q[k,])) {
        w <- Z[i,i] - Z[i,k] # Money pump cost
        b <- Z[i,i]          # Budget
        g_min[[i]][[length(g_min[[i]]) + 1]] <- c(k, w, b)
      }
    }
  }
  
  # Graph for Maximum MPI (flips weights)
  g_max <- vector("list", t)
  for(i in 1:t) {
    for(k in 1:t) {
      if(k != i && G[i,k] == 1 && !all(Q[i,] == Q[k,])) {
        w <- Z[i,k] - Z[i,i]
        b <- Z[i,i]
        g_max[[i]][[length(g_max[[i]]) + 1]] <- c(k, w, b)
      }
    }
  }
  
  # If perfectly rational, bounds are zero
  if(cycle_detection_topo(g_min) == 0) {
    return(c(0, 0))
  }
  
  # Binary search to find exact minimum and maximum MPIs
  min_MPI <- minimum_cost_time(g_min, -10000, 10000)
  max_MPI <- -1 * minimum_cost_time(g_max, -10000, 10000)
  return(c(min_MPI, max_MPI))
}
