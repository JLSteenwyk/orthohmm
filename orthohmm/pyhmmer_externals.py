import os
import sys
from typing import List, Dict, Tuple
import numpy as np

import pyhmmer
from pyhmmer import easel

from .helpers import SubstitutionMatrix


def execute_pyhmmer_search_batch_optimized(
    files: List[str],
    fasta_directory: str,
    substitution_matrix: SubstitutionMatrix,
    evalue_threshold: float,
    cpu: int,
) -> Dict[Tuple[str, str], List[Dict]]:
    """Optimized PyHMMER implementation using batch processing"""
    
    import itertools
    
    alphabet = easel.Alphabet.amino()
    
    # Load all sequences once and keep them in memory
    sequences_cache = {}
    file_to_seq_names = {}
    
    print("   Loading all sequences into memory...")
    for file in files:
        file_path = os.path.join(fasta_directory, file)
        with easel.SequenceFile(file_path, digital=True, alphabet=alphabet) as seq_file:
            seqs = seq_file.read_block()
            sequences_cache[file] = seqs
            file_to_seq_names[file] = [seq.name.decode() for seq in seqs]
    
    all_results = {}
    
    # Process each file pair
    file_pairs = list(itertools.product(files, repeat=2))
    print(f"   Running {len(file_pairs)} pairwise searches with batch processing...")
    
    for i, (query_file, target_file) in enumerate(file_pairs):
        progress = ((i + 1) / len(file_pairs)) * 100
        sys.stdout.write(f"\r          {progress:.1f}% complete")
        sys.stdout.flush()
        
        query_sequences = sequences_cache[query_file]
        target_sequences = sequences_cache[target_file] 
        query_names = file_to_seq_names[query_file]
        
        results = []
        
        # Use pyhmmer's efficient batch processing
        try:
            hit_iter = pyhmmer.phmmer(query_sequences, target_sequences, cpus=cpu)
            
            for query_idx, hits in enumerate(hit_iter):
                query_name = query_names[query_idx]
                
                # Process hits for this query
                for hit_idx in range(len(hits)):
                    hit = hits[hit_idx]
                    if hit.evalue <= evalue_threshold:
                        results.append({
                            'target_name': hit.name.decode(),
                            'query_name': query_name,
                            'evalue': hit.evalue,
                            'score': hit.score
                        })
        except Exception as e:
            print(f"\nError processing {query_file} vs {target_file}: {e}")
            results = []
        
        all_results[(query_file, target_file)] = results
    
    print()  # New line after progress
    return all_results


def pyhmmer_results_to_numpy(results: List[Dict]) -> np.ndarray:
    """Convert pyhmmer results to numpy array format expected by helpers.py"""
    if not results:
        # Return empty array with correct dtype
        dtype_res = [
            ("target_name", "U50"),
            ("query_name", "U50"), 
            ("evalue", float),
            ("score", float)
        ]
        return np.array([], dtype=dtype_res)
    
    # Convert to numpy structured array
    dtype_res = [
        ("target_name", "U50"),
        ("query_name", "U50"),
        ("evalue", float),
        ("score", float)
    ]
    
    numpy_data = np.empty(len(results), dtype=dtype_res)
    for i, result in enumerate(results):
        numpy_data[i] = (
            result['target_name'],
            result['query_name'],
            result['evalue'],
            result['score']
        )
    
    return numpy_data