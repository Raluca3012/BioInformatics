#Formulate 3 scoring ecuations that are able to show the level of similarity
#between the 2 seq aligned in the previous assignment
#Implement each of these scoring ecuations in your current implementations

import os
import math
import time
import matplotlib.pyplot as plt
import numpy as np

PATH_INFLUENZA = "influenza_ref.fasta" 
PATH_COVID = "covid_ref.fasta"
MATCH_REWARD = 3
MISMATCH_PENALTY = -3
GAP_PENALTY = -2

def read_fasta_native(filepath):
    if not os.path.exists(filepath): return None
    sequence = ""
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith(">"): sequence += line.strip().upper()
    return sequence

def smith_waterman_kernel_raw(seq1, seq2):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0] * cols for _ in range(rows)]
    max_score = 0
    
    for i in range(1, rows):
        for j in range(1, cols):
            match_op = matrix[i-1][j-1] + (MATCH_REWARD if seq1[i-1] == seq2[j-1] else MISMATCH_PENALTY)
            delete_op = matrix[i-1][j] + GAP_PENALTY
            insert_op = matrix[i][j-1] + GAP_PENALTY

            cell_score = max(0, match_op, delete_op, insert_op)
            matrix[i][j] = cell_score
            
            if cell_score > max_score:
                max_score = cell_score
                
    return max_score

def run_simulation_with_metrics(genome_a, genome_b, window_size=150, step=100):
    len_a, len_b = len(genome_a), len(genome_b)

    windows_a = [genome_a[i:i+window_size] for i in range(0, len_a, step) if len(genome_a[i:i+window_size]) > window_size//2]
    windows_b = [genome_b[i:i+window_size] for i in range(0, len_b, step) if len(genome_b[i:i+window_size]) > window_size//2]
    
    print(f"Windows: A={len(windows_a)}, B={len(windows_b)}")

    raw_scores_flat = []
    matrix_raw = np.zeros((len(windows_a), len(windows_b)))
    
    print("Calculating Raw Scores (Equation 1)...")
    for i, wa in enumerate(windows_a):
        for j, wb in enumerate(windows_b):
            score = smith_waterman_kernel_raw(wa, wb)
            matrix_raw[i][j] = score
            raw_scores_flat.append(score)
            
    mean_score = sum(raw_scores_flat) / len(raw_scores_flat)
    variance = sum([(x - mean_score) ** 2 for x in raw_scores_flat]) / len(raw_scores_flat)
    std_dev = math.sqrt(variance)

    matrix_norm = np.zeros(matrix_raw.shape)
    matrix_z = np.zeros(matrix_raw.shape)
    
    print("Calculating Normalized & Z-Scores (Equations 2 & 3)...")
    max_possible_score = window_size * MATCH_REWARD
    
    for i in range(len(windows_a)):
        for j in range(len(windows_b)):
            raw_s = matrix_raw[i][j]

            pct_score = (raw_s / max_possible_score) * 100
            matrix_norm[i][j] = pct_score
   
            if std_dev > 0:
                z_score = (raw_s - mean_score) / std_dev
            else:
                z_score = 0
            matrix_z[i][j] = z_score

    return matrix_raw, matrix_norm, matrix_z

if __name__ == "__main__":

    seq_flu = read_fasta_native(PATH_INFLUENZA)
    seq_cov = read_fasta_native(PATH_COVID)

    mat_raw, mat_norm, mat_z = run_simulation_with_metrics(seq_flu, seq_cov, window_size=150, step=100)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    im1 = axes[0].imshow(mat_raw, cmap='inferno', aspect='auto')
    axes[0].set_title("Eq 1: Raw Alignment Scores (S)\n(Absolute similarity)")
    plt.colorbar(im1, ax=axes[0], label="Score")

    im2 = axes[1].imshow(mat_norm, cmap='viridis', aspect='auto')
    axes[1].set_title("Eq 2: Normalized Similarity (%)\n(Relative to window size)")
    plt.colorbar(im2, ax=axes[1], label="Similarity %")
    
    im3 = axes[2].imshow(mat_z, cmap='coolwarm', aspect='auto')
    axes[2].set_title("Eq 3: Z-Scores (Significance)\n(Deviations from mean)")
    plt.colorbar(im3, ax=axes[2], label="Standard Deviations")
    
    plt.tight_layout()
    plt.show()
    
    print("Implementation complete.")
    print(f"Max Raw Score: {np.max(mat_raw)}")
    print(f"Max Similarity: {np.max(mat_norm):.2f}%")
    print(f"Max Z-Score: {np.max(mat_z):.2f}")