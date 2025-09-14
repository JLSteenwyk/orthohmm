"""
Alternative clustering methods for OrthoHMM
Provides different approaches to cluster protein similarity networks into orthogroups
"""

import numpy as np
import subprocess
import tempfile
import os
from enum import Enum
from typing import List, Dict, Tuple, Optional
import warnings

try:
    import igraph as ig
    import leidenalg
    LEIDEN_AVAILABLE = True
except ImportError:
    LEIDEN_AVAILABLE = False
    ig = None
    leidenalg = None

try:
    from sklearn.cluster import SpectralClustering, AffinityPropagation
    from sklearn.metrics.pairwise import rbf_kernel
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    SpectralClustering = None
    AffinityPropagation = None


class ClusteringMethod(Enum):
    """Available clustering methods"""
    MCL = "mcl"                        # Markov Cluster Algorithm (default)
    LEIDEN = "leiden"                  # Leiden algorithm for community detection
    SPECTRAL = "spectral"              # Spectral clustering
    AFFINITY_PROPAGATION = "affinity_propagation"  # Affinity Propagation


def cluster_network(
    edges_file: str,
    output_file: str,
    method: ClusteringMethod = ClusteringMethod.MCL,
    mcl_path: str = "mcl",
    inflation: float = 1.5,
    cpu: int = 1,
    **kwargs
) -> None:
    """
    Cluster protein similarity network using specified method
    
    Args:
        edges_file: Path to edges file (tab-separated: protein1, protein2, score)
        output_file: Path to output clustered results
        method: Clustering method to use
        mcl_path: Path to MCL binary (for MCL method)
        inflation: MCL inflation parameter (for MCL method)
        cpu: Number of CPU cores to use
        **kwargs: Additional method-specific parameters
    """
    
    if method == ClusteringMethod.MCL:
        _cluster_with_mcl(edges_file, output_file, mcl_path, inflation, cpu)
    elif method == ClusteringMethod.LEIDEN:
        if not LEIDEN_AVAILABLE:
            raise ImportError("Leiden clustering requires 'igraph' and 'leidenalg' packages. Install with: pip install igraph leidenalg")
        _cluster_with_leiden(edges_file, output_file, **kwargs)
    elif method == ClusteringMethod.SPECTRAL:
        if not SKLEARN_AVAILABLE:
            raise ImportError("Spectral clustering requires 'scikit-learn' package. Install with: pip install scikit-learn")
        _cluster_with_spectral(edges_file, output_file, **kwargs)
    elif method == ClusteringMethod.AFFINITY_PROPAGATION:
        if not SKLEARN_AVAILABLE:
            raise ImportError("Affinity Propagation requires 'scikit-learn' package. Install with: pip install scikit-learn")
        _cluster_with_affinity_propagation(edges_file, output_file, **kwargs)
    else:
        raise ValueError(f"Unknown clustering method: {method}")


def _cluster_with_mcl(
    edges_file: str,
    output_file: str,
    mcl_path: str,
    inflation: float,
    cpu: int
) -> None:
    """Original MCL clustering implementation"""
    cmd = f"{mcl_path} {edges_file} -te {cpu} --abc -I {inflation} -o {output_file}"
    subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _cluster_with_leiden(
    edges_file: str,
    output_file: str,
    resolution: float = 1.0,
    random_state: Optional[int] = None,
    **kwargs
) -> None:
    """Cluster using Leiden algorithm"""
    
    # Read edges and build graph
    edges, weights, node_map = _read_edges_file(edges_file)
    
    # Create igraph Graph
    g = ig.Graph()
    g.add_vertices(len(node_map))
    g.add_edges(edges)
    g.es['weight'] = weights
    
    # Set vertex names
    vertex_names = [''] * len(node_map)
    for name, idx in node_map.items():
        vertex_names[idx] = name
    g.vs['name'] = vertex_names
    
    # Run Leiden clustering
    if random_state is not None:
        ig.set_random_number_generator(ig.RNG_GSL)
        ig.Graph.rng.set_seed(random_state)
    
    partition = leidenalg.find_partition(
        g, 
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        weights='weight'
    )
    
    # Write results in MCL format
    _write_clusters_mcl_format(partition, vertex_names, output_file)


def _cluster_with_spectral(
    edges_file: str,
    output_file: str,
    n_clusters: Optional[int] = None,
    gamma: float = 1.0,
    random_state: Optional[int] = None,
    **kwargs
) -> None:
    """Cluster using Spectral clustering"""
    
    # Read edges and build similarity matrix
    edges, weights, node_map = _read_edges_file(edges_file)
    n_nodes = len(node_map)
    
    # Build adjacency matrix
    adj_matrix = np.zeros((n_nodes, n_nodes))
    for (i, j), weight in zip(edges, weights):
        adj_matrix[i, j] = weight
        adj_matrix[j, i] = weight  # Make symmetric
    
    # Estimate number of clusters if not provided
    if n_clusters is None:
        n_clusters = _estimate_n_clusters(adj_matrix)
    
    # Apply spectral clustering
    spectral = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed',
        random_state=random_state,
        n_jobs=-1
    )
    
    cluster_labels = spectral.fit_predict(adj_matrix)
    
    # Convert to MCL format
    node_names = [''] * len(node_map)
    for name, idx in node_map.items():
        node_names[idx] = name
    
    _write_clusters_from_labels(cluster_labels, node_names, output_file)


def _cluster_with_affinity_propagation(
    edges_file: str,
    output_file: str,
    damping: float = 0.5,
    preference: Optional[float] = None,
    random_state: Optional[int] = None,
    **kwargs
) -> None:
    """Cluster using Affinity Propagation"""
    
    # Read edges and build similarity matrix
    edges, weights, node_map = _read_edges_file(edges_file)
    n_nodes = len(node_map)
    
    # Build similarity matrix (negative distances for AP)
    similarity_matrix = np.full((n_nodes, n_nodes), -np.inf)
    np.fill_diagonal(similarity_matrix, 0)
    
    for (i, j), weight in zip(edges, weights):
        similarity_matrix[i, j] = weight
        similarity_matrix[j, i] = weight
    
    # Set preference (median similarity if not provided)
    if preference is None:
        finite_similarities = similarity_matrix[similarity_matrix != -np.inf]
        if len(finite_similarities) > 0:
            preference = np.median(finite_similarities)
        else:
            preference = 0.0
    
    # Apply Affinity Propagation
    ap = AffinityPropagation(
        affinity='precomputed',
        damping=damping,
        preference=preference,
        random_state=random_state
    )
    
    cluster_labels = ap.fit_predict(similarity_matrix)
    
    # Convert to MCL format
    node_names = [''] * len(node_map)
    for name, idx in node_map.items():
        node_names[idx] = name
    
    _write_clusters_from_labels(cluster_labels, node_names, output_file)


def _read_edges_file(edges_file: str) -> Tuple[List[Tuple[int, int]], List[float], Dict[str, int]]:
    """Read edges file and return edges, weights, and node mapping"""
    
    edges = []
    weights = []
    node_map = {}
    node_counter = 0
    
    with open(edges_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    node1, node2, weight = parts[0], parts[1], float(parts[2])
                    
                    # Map nodes to integers
                    if node1 not in node_map:
                        node_map[node1] = node_counter
                        node_counter += 1
                    if node2 not in node_map:
                        node_map[node2] = node_counter
                        node_counter += 1
                    
                    edges.append((node_map[node1], node_map[node2]))
                    weights.append(weight)
    
    return edges, weights, node_map


def _estimate_n_clusters(adj_matrix: np.ndarray) -> int:
    """Estimate number of clusters using eigenvalue gap heuristic"""
    
    # Compute Laplacian eigenvalues
    try:
        from scipy.linalg import eigvals
        degree = np.sum(adj_matrix, axis=1)
        laplacian = np.diag(degree) - adj_matrix
        eigenvals = eigvals(laplacian)
        eigenvals = np.sort(eigenvals.real)
        
        # Find largest eigenvalue gap
        gaps = np.diff(eigenvals)
        n_clusters = np.argmax(gaps) + 2  # +2 because we want clusters, not gaps
        
        # Reasonable bounds
        n_clusters = max(2, min(n_clusters, len(adj_matrix) // 10))
        
    except ImportError:
        # Fallback: estimate based on graph density
        n_edges = np.sum(adj_matrix > 0) // 2
        n_nodes = adj_matrix.shape[0]
        density = n_edges / (n_nodes * (n_nodes - 1) / 2) if n_nodes > 1 else 0
        n_clusters = max(2, min(int(n_nodes * density * 10), n_nodes // 5))
    
    return n_clusters


def _write_clusters_mcl_format(partition, vertex_names: List[str], output_file: str) -> None:
    """Write igraph partition in MCL output format"""
    
    with open(output_file, 'w') as f:
        for cluster in partition:
            cluster_proteins = [vertex_names[i] for i in cluster]
            if cluster_proteins:  # Skip empty clusters
                f.write('\t'.join(cluster_proteins) + '\n')


def _write_clusters_from_labels(cluster_labels: np.ndarray, node_names: List[str], output_file: str) -> None:
    """Write clustering results from cluster labels in MCL output format"""
    
    # Group nodes by cluster label
    clusters = {}
    for i, label in enumerate(cluster_labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(node_names[i])
    
    # Write results
    with open(output_file, 'w') as f:
        for cluster_id, proteins in clusters.items():
            if len(proteins) > 0:  # Skip empty clusters
                f.write('\t'.join(proteins) + '\n')


def get_method_description(method: ClusteringMethod) -> str:
    """Get human-readable description of clustering method"""
    descriptions = {
        ClusteringMethod.MCL: "Markov Cluster Algorithm - Flow-based graph clustering",
        ClusteringMethod.LEIDEN: "Leiden Algorithm - Community detection with modularity optimization", 
        ClusteringMethod.SPECTRAL: "Spectral Clustering - Eigendecomposition-based clustering",
        ClusteringMethod.AFFINITY_PROPAGATION: "Affinity Propagation - Exemplar-based clustering"
    }
    return descriptions.get(method, "Unknown method")


def demonstrate_clustering_methods():
    """Demonstrate different clustering methods with example data"""
    
    print("=== Clustering Method Comparison ===\n")
    
    # Create example edges file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_edges:
        # Simple example network
        edges_data = [
            "protein1\tprotein2\t0.9",
            "protein1\tprotein3\t0.8", 
            "protein2\tprotein3\t0.85",
            "protein4\tprotein5\t0.95",
            "protein4\tprotein6\t0.7",
            "protein5\tprotein6\t0.8",
            "protein7\tprotein8\t0.6"
        ]
        tmp_edges.write('\n'.join(edges_data))
        tmp_edges_path = tmp_edges.name
    
    print("Example Network:")
    for edge in edges_data:
        print(f"  {edge}")
    print()
    
    # Test available methods
    available_methods = [ClusteringMethod.MCL]
    if LEIDEN_AVAILABLE:
        available_methods.append(ClusteringMethod.LEIDEN)
    if SKLEARN_AVAILABLE:
        available_methods.extend([ClusteringMethod.SPECTRAL, ClusteringMethod.AFFINITY_PROPAGATION])
    
    print("Clustering Results:")
    for method in available_methods:
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_out:
                tmp_out_path = tmp_out.name
            
            print(f"\n{method.value.upper()}:")
            print(f"  Description: {get_method_description(method)}")
            
            if method == ClusteringMethod.MCL:
                print("  Note: Requires MCL binary to be installed")
                continue
            else:
                cluster_network(tmp_edges_path, tmp_out_path, method)
                
                with open(tmp_out_path, 'r') as f:
                    clusters = f.readlines()
                    for i, cluster in enumerate(clusters):
                        proteins = cluster.strip().split('\t')
                        print(f"    Cluster {i+1}: {', '.join(proteins)}")
                
                os.unlink(tmp_out_path)
                
        except Exception as e:
            print(f"    Error: {e}")
    
    # Cleanup
    os.unlink(tmp_edges_path)


if __name__ == "__main__":
    demonstrate_clustering_methods()