#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <queue>
#include <algorithm>

using namespace Rcpp;

struct Edge {
    int to;
    double w;
    double b;
};

typedef std::vector<std::vector<Edge>> Graph;

// cycle_detection_topo
int cycle_detection_topo(const Graph& g) {
    int n = g.size();
    std::vector<int> indegree(n, 0);
    
    for (int i = 0; i < n; ++i) {
        for (const auto& edge : g[i]) {
            indegree[edge.to]++;
        }
    }
    
    std::queue<int> q;
    for (int i = 0; i < n; ++i) {
        if (indegree[i] == 0) {
            q.push(i);
        }
    }
    
    int count = 0;
    while (!q.empty()) {
        int node = q.front();
        q.pop();
        count++;
        
        for (const auto& edge : g[node]) {
            indegree[edge.to]--;
            if (indegree[edge.to] == 0) {
                q.push(edge.to);
            }
        }
    }
    
    if (count < n) return 1; // cycle exists
    return 0; // acyclic
}

// bellmanford
// Returns false if negative cycle found, true if successful
bool bellmanford(const Graph& g, int s, std::vector<double>& d) {
    // Initialize with virtual supersource distance of 0.0 to all vertices
    int n = g.size();
    d.assign(n, 0.0);
    
    for (int t = 0; t < n - 1; ++t) {
        for (int i = 0; i < n; ++i) {
            for (const auto& edge : g[i]) {
                if (d[i] + edge.w < d[edge.to]) {
                    d[edge.to] = d[i] + edge.w;
                }
            }
        }
    }
    
    for (int i = 0; i < n; ++i) {
        for (const auto& edge : g[i]) {
            if (d[i] + edge.w < d[edge.to]) {
                return false; // Negative cycle
            }
        }
    }
    
    return true;
}

// zero_residual_graph
Graph zero_residual_graph(const Graph& g, const std::vector<double>& d) {
    int n = g.size();
    Graph e(n);
    for (int i = 0; i < n; ++i) {
        for (const auto& edge : g[i]) {
            if (std::abs(edge.w + d[i] - d[edge.to]) < 1e-9) {
                e[i].push_back({edge.to, 0.0, 0.0});
            }
        }
    }
    return e;
}

// redefined_graph
Graph redefined_graph(const Graph& g_init, double u) {
    Graph g = g_init;
    int n = g.size();
    for (int i = 0; i < n; ++i) {
        for (auto& edge : g[i]) {
            edge.w = edge.w - u * edge.b;
        }
    }
    return g;
}

// minimum_cost_time
double minimum_cost_time(const Graph& g_init, double low, double high) {
    if (cycle_detection_topo(g_init) == 0) return 0.0;
    
    while (std::abs(high - low) > 1e-10) {
        double u = (low + high) / 2.0;
        Graph g = redefined_graph(g_init, u);
        
        std::vector<double> d;
        bool success = bellmanford(g, 0, d);
        
        if (success) {
            Graph e = zero_residual_graph(g, d);
            if (cycle_detection_topo(e) == 1) {
                return u;
            } else {
                low = u;
            }
        } else {
            high = u;
        }
    }
    return (high + low) / 2.0;
}

// [[Rcpp::export]]
NumericVector compute_min_max_mpi_cpp(NumericMatrix Z, NumericMatrix Q) {
    int t = Z.nrow();
    
    Graph g_min(t);
    Graph g_max(t);
    
    for (int i = 0; i < t; ++i) {
        for (int k = 0; k < t; ++k) {
            if (k != i && Z(i, k) <= Z(i, i)) {
                // Check if all(Q[i,] == Q[k,])
                bool same_q = true;
                for (int c = 0; c < Q.ncol(); ++c) {
                    if (Q(i, c) != Q(k, c)) {
                        same_q = false;
                        break;
                    }
                }
                
                if (!same_q) {
                    double w_min = Z(i, i) - Z(i, k);
                    double b_min = Z(i, i);
                    g_min[i].push_back({k, w_min, b_min});
                    
                    double w_max = Z(i, k) - Z(i, i);
                    double b_max = Z(i, i);
                    g_max[i].push_back({k, w_max, b_max});
                }
            }
        }
    }
    
    if (cycle_detection_topo(g_min) == 0) {
        return NumericVector::create(0.0, 0.0);
    }
    
    double min_MPI = minimum_cost_time(g_min, -10000.0, 10000.0);
    double max_MPI = -1.0 * minimum_cost_time(g_max, -10000.0, 10000.0);
    
    return NumericVector::create(min_MPI, max_MPI);
}
