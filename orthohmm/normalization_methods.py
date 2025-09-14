"""
Alternative gene length normalization methods for OrthoHMM
Provides different approaches to account for protein length bias in HMMER scores
"""

import numpy as np
import math
from enum import Enum
from typing import Union


class NormalizationMethod(Enum):
    """Available normalization methods"""
    DEFAULT = "default"                    # raw_score / (target_len + query_len)
    GEOMETRIC_MEAN = "geometric_mean"      # raw_score / sqrt(target_len * query_len)  
    MAX_LENGTH = "max_length"              # raw_score / max(target_len, query_len)
    LOGARITHMIC = "logarithmic"            # raw_score / log(target_len + query_len + 1)
    MIN_LENGTH = "min_length"              # raw_score / min(target_len, query_len)


def normalize_by_gene_length_method(
    res_merged: np.ndarray, 
    method: NormalizationMethod = NormalizationMethod.DEFAULT
) -> np.ndarray:
    """
    Normalize scores by gene length using specified method
    
    Args:
        res_merged: Array with columns [target_name, query_name, evalue, score, target_len, query_len]
        method: Normalization method to use
        
    Returns:
        Modified array with normalized scores in column 3
    """
    
    # Extract lengths and scores
    scores = res_merged[:, 3].astype(float)
    target_lengths = res_merged[:, 4].astype(float)
    query_lengths = res_merged[:, 5].astype(float)
    
    if method == NormalizationMethod.DEFAULT:
        # Current default: arithmetic mean
        normalized_scores = scores / (target_lengths + query_lengths)
        
    elif method == NormalizationMethod.GEOMETRIC_MEAN:
        # Geometric mean normalization
        geometric_means = np.sqrt(target_lengths * query_lengths)
        # Avoid division by zero
        geometric_means = np.where(geometric_means == 0, 1, geometric_means)
        normalized_scores = scores / geometric_means
        
    elif method == NormalizationMethod.MAX_LENGTH:
        # Normalize by longer sequence
        max_lengths = np.maximum(target_lengths, query_lengths)
        # Avoid division by zero
        max_lengths = np.where(max_lengths == 0, 1, max_lengths)
        normalized_scores = scores / max_lengths
        
    elif method == NormalizationMethod.LOGARITHMIC:
        # Logarithmic scaling
        log_lengths = np.log(target_lengths + query_lengths + 1)  # +1 prevents log(0)
        normalized_scores = scores / log_lengths
        
    elif method == NormalizationMethod.MIN_LENGTH:
        # Normalize by shorter sequence
        min_lengths = np.minimum(target_lengths, query_lengths)
        # Avoid division by zero
        min_lengths = np.where(min_lengths == 0, 1, min_lengths)
        normalized_scores = scores / min_lengths
        
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    # Update the array
    res_merged[:, 3] = normalized_scores
    return res_merged


def compare_normalization_methods(
    res_merged: np.ndarray,
    sample_size: int = 5
) -> dict:
    """
    Compare different normalization methods on sample data
    
    Args:
        res_merged: Array with HMMER results and gene lengths
        sample_size: Number of examples to show
        
    Returns:
        Dictionary with comparison results
    """
    
    if len(res_merged) == 0:
        return {"error": "No data provided"}
    
    # Take a sample for comparison
    sample_indices = np.random.choice(len(res_merged), min(sample_size, len(res_merged)), replace=False)
    sample_data = res_merged[sample_indices].copy()
    
    results = {
        "sample_size": len(sample_data),
        "methods": {},
        "original_data": []
    }
    
    # Store original data info
    for i, row in enumerate(sample_data):
        target, query, evalue, score, tlen, qlen = row
        results["original_data"].append({
            "index": i,
            "query": query,
            "target": target,
            "raw_score": float(score),
            "target_length": int(tlen),
            "query_length": int(qlen),
            "combined_length": int(tlen) + int(qlen)
        })
    
    # Test each normalization method
    for method in NormalizationMethod:
        test_data = sample_data.copy()
        normalized_data = normalize_by_gene_length_method(test_data, method)
        
        results["methods"][method.value] = {
            "name": method.value,
            "description": get_method_description(method),
            "normalized_scores": [float(score) for score in normalized_data[:, 3]]
        }
    
    return results


def get_method_description(method: NormalizationMethod) -> str:
    """Get human-readable description of normalization method"""
    descriptions = {
        NormalizationMethod.DEFAULT: "score / (target_length + query_length)",
        NormalizationMethod.GEOMETRIC_MEAN: "score / sqrt(target_length × query_length)", 
        NormalizationMethod.MAX_LENGTH: "score / max(target_length, query_length)",
        NormalizationMethod.LOGARITHMIC: "score / log(target_length + query_length + 1)",
        NormalizationMethod.MIN_LENGTH: "score / min(target_length, query_length)"
    }
    return descriptions.get(method, "Unknown method")


def demonstrate_normalization_methods():
    """Demonstrate different normalization methods with example data"""
    
    print("=== Gene Length Normalization Method Comparison ===\\n")
    
    # Create example data with different protein length scenarios
    example_data = np.array([
        ['short_A', 'short_B', 1e-50, 150.0, 200, 250],      # Short proteins
        ['medium_A', 'medium_B', 1e-45, 300.0, 500, 600],    # Medium proteins
        ['long_A', 'long_B', 1e-40, 450.0, 1000, 1200],      # Long proteins
        ['mixed_A', 'mixed_B', 1e-35, 250.0, 300, 1000],     # Very different lengths
    ], dtype=object)
    
    # Convert numeric columns
    example_data[:, 2] = example_data[:, 2].astype(float)  # evalue
    example_data[:, 3] = example_data[:, 3].astype(float)  # score
    example_data[:, 4] = example_data[:, 4].astype(int)    # target_length
    example_data[:, 5] = example_data[:, 5].astype(int)    # query_length
    
    # Show original data
    print("Original Data:")
    for i, row in enumerate(example_data):
        target, query, evalue, score, tlen, qlen = row
        print(f"  {i+1}. {query} → {target}: Score={score}, Lengths=({tlen}, {qlen})")
    print()
    
    # Compare methods
    print("Normalization Results:")
    print(f"{'Method':<20} {'Short':<12} {'Medium':<12} {'Long':<12} {'Mixed':<12}")
    print("-" * 70)
    
    for method in NormalizationMethod:
        test_data = example_data.copy()
        normalized = normalize_by_gene_length_method(test_data, method)
        scores = [f"{score:.6f}" for score in normalized[:, 3].astype(float)]
        
        print(f"{method.value:<20} {scores[0]:<12} {scores[1]:<12} {scores[2]:<12} {scores[3]:<12}")
    
    print()
    print("Method Descriptions:")
    for method in NormalizationMethod:
        print(f"  {method.value}: {get_method_description(method)}")


if __name__ == "__main__":
    demonstrate_normalization_methods()