"""
GPU acceleration utilities for OrthoHMM using CuPy
Provides GPU-accelerated matrix operations with CPU fallback
"""

import numpy as np
from typing import Any, Union

# Try to import CuPy, fallback to NumPy if not available
try:
    import cupy as cp
    GPU_AVAILABLE = True
    
    # Test if GPU is actually accessible
    try:
        test_array = cp.array([1, 2, 3])
        _ = cp.asnumpy(test_array)  # Force GPU operation
        GPU_FUNCTIONAL = True
    except Exception:
        GPU_FUNCTIONAL = False
        GPU_AVAILABLE = False
        
except ImportError:
    # CuPy not installed, use NumPy as fallback
    import numpy as cp
    GPU_AVAILABLE = False
    GPU_FUNCTIONAL = False


def get_array_module(array: Any) -> Any:
    """
    Get the appropriate array module (cupy or numpy) for the given array
    
    Args:
        array: Input array (numpy or cupy)
        
    Returns:
        cp (CuPy) if GPU available and array is on GPU, otherwise np (NumPy)
    """
    if GPU_AVAILABLE and hasattr(cp, 'get_array_module'):
        return cp.get_array_module(array)
    else:
        return np


def to_gpu(array: np.ndarray) -> Union[np.ndarray, Any]:
    """
    Move array to GPU if available, otherwise return original array
    
    Args:
        array: NumPy array
        
    Returns:
        CuPy array if GPU available, otherwise original NumPy array
    """
    if GPU_FUNCTIONAL:
        try:
            return cp.asarray(array)
        except Exception:
            return array
    return array


def to_cpu(array: Union[np.ndarray, Any]) -> np.ndarray:
    """
    Move array to CPU (convert CuPy array to NumPy)
    
    Args:
        array: NumPy or CuPy array
        
    Returns:
        NumPy array
    """
    if GPU_AVAILABLE and hasattr(cp, 'asnumpy'):
        try:
            return cp.asnumpy(array)
        except Exception:
            return np.asarray(array)
    return np.asarray(array)


def gpu_info() -> dict:
    """
    Get GPU availability and status information
    
    Returns:
        Dictionary with GPU status information
    """
    info = {
        'cupy_installed': GPU_AVAILABLE,
        'gpu_functional': GPU_FUNCTIONAL,
        'device_count': 0,
        'device_name': None,
        'memory_info': None
    }
    
    if GPU_FUNCTIONAL:
        try:
            info['device_count'] = cp.cuda.runtime.getDeviceCount()
            device = cp.cuda.Device()
            info['device_name'] = device.attributes['Name']
            mem_info = cp.cuda.MemoryInfo()
            info['memory_info'] = {
                'total': mem_info.total,
                'free': mem_info.free,
                'used': mem_info.total - mem_info.free
            }
        except Exception as e:
            info['error'] = str(e)
    
    return info


def normalize_by_gene_length_gpu(res_merged: np.ndarray) -> np.ndarray:
    """
    GPU-accelerated version of normalize_by_gene_length
    
    Args:
        res_merged: Array with columns [target_name, query_name, evalue, score, target_len, query_len]
        
    Returns:
        Modified array with normalized scores
    """
    if GPU_FUNCTIONAL and len(res_merged) > 100:  # Only use GPU for larger arrays
        try:
            # Move to GPU
            gpu_array = cp.asarray(res_merged)
            
            # Perform normalization on GPU
            # res_merged[:, 3] = res_merged[:, 3] / (res_merged[:, 4] + res_merged[:, 5])
            gpu_array[:, 3] = gpu_array[:, 3] / (gpu_array[:, 4] + gpu_array[:, 5])
            
            # Move back to CPU
            return cp.asnumpy(gpu_array)
            
        except Exception:
            # Fallback to CPU if GPU operation fails
            pass
    
    # CPU fallback (original implementation)
    res_merged[:, 3] = res_merged[:, 3] / (res_merged[:, 4] + res_merged[:, 5])
    return res_merged


def batch_score_calculations_gpu(
    scores: np.ndarray, 
    query_lengths: np.ndarray, 
    target_lengths: np.ndarray,
    correction_factors: np.ndarray
) -> np.ndarray:
    """
    GPU-accelerated batch scoring calculations
    
    Args:
        scores: Array of raw scores
        query_lengths: Array of query sequence lengths
        target_lengths: Array of target sequence lengths  
        correction_factors: Array of phylogenetic correction factors
        
    Returns:
        Array of normalized scores
    """
    if GPU_FUNCTIONAL and len(scores) > 1000:  # Use GPU for large batches
        try:
            # Move all arrays to GPU
            gpu_scores = cp.asarray(scores)
            gpu_query_len = cp.asarray(query_lengths)
            gpu_target_len = cp.asarray(target_lengths)
            gpu_corrections = cp.asarray(correction_factors)
            
            # Vectorized calculation on GPU
            normalized_scores = (gpu_scores / (gpu_query_len + gpu_target_len)) / gpu_corrections
            
            # Move result back to CPU
            return cp.asnumpy(normalized_scores)
            
        except Exception:
            # Fallback to CPU
            pass
    
    # CPU fallback
    return (scores / (query_lengths + target_lengths)) / correction_factors


def print_gpu_status():
    """Print GPU status information"""
    info = gpu_info()
    
    if info['cupy_installed']:
        if info['gpu_functional']:
            print(f"🚀 GPU acceleration: ENABLED")
            if info['device_name']:
                print(f"   Device: {info['device_name']}")
            if info['memory_info']:
                mem = info['memory_info']
                print(f"   Memory: {mem['used']/1e9:.1f}GB / {mem['total']/1e9:.1f}GB used")
        else:
            print("⚠️  GPU acceleration: DISABLED (CuPy installed but GPU not functional)")
    else:
        print("💻 GPU acceleration: DISABLED (CuPy not installed, using CPU)")