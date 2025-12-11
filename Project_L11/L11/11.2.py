#Download from NCBI the influenza genome and covid-19 genomes and align these
#2 genomes by using the local alignment method
#NOTE that you have to add in between layers solutions in order to be able to 
#align the 2 genome files
#HINT The alignment order setp-by-step on big regions with the connection
#between the results
#Note that the sequence align algorithm does not allow for align of sequences
#because the bigger seq. the bigger the scoring matrix
#The main result should be the vizualisation of the simulations between the 2
#genomes
#NNOTE Please do not use the shortcuts given by AI, use NATIVE CODE 

import os
import math
import time
import matplotlib.pyplot as plt
import numpy as np 
PATH_INFLUENZA = "influenza_ref.fasta" 
PATH_COVID = "covid_ref.fasta"

def read_fasta_native(filepath):
    sequence = ""
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line.upper()
    return sequence

def smith_waterman_kernel(seq1, seq2, match=3, mismatch=-3, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    max_score = 0
    
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i-1] == seq2[j-1]:
                score_diag = matrix[i-1][j-1] + match
            else:
                score_diag = matrix[i-1][j-1] + mismatch
          
            score_up = matrix[i-1][j] + gap
            score_left = matrix[i][j-1] + gap
            current_score = max(0, score_diag, score_up, score_left)
            matrix[i][j] = current_score
            if current_score > max_score:
                max_score = current_score
                
    return max_score

def layered_alignment_simulation(genome_a, genome_b, window_size=100, step=50):
    len_a = len(genome_a)
    len_b = len(genome_b)
    
    print(f"Length Influenza: {len_a} bp")
    print(f"Length Covid: {len_b} bp")

    windows_a = [genome_a[i:i+window_size] for i in range(0, len_a, step) if len(genome_a[i:i+window_size]) > window_size//2]
    windows_b = [genome_b[i:i+window_size] for i in range(0, len_b, step) if len(genome_b[i:i+window_size]) > window_size//2]

    similarity_map = np.zeros((len(windows_a), len(windows_b)))
    start_time = time.time()
    
    for i, seq_a in enumerate(windows_a):
        for j, seq_b in enumerate(windows_b):
            score = smith_waterman_kernel(seq_a, seq_b)
            similarity_map[i][j] = score
            
        if i % 5 == 0:
            print(f"Processed row {i}/{len(windows_a)}...")

    print(f"Done!")
    return similarity_map

if __name__ == "__main__":
    seq_flu = read_fasta_native(PATH_INFLUENZA)   
    seq_cov = read_fasta_native(PATH_COVID)
    result_matrix = layered_alignment_simulation(seq_flu, seq_cov, window_size=150, step=100)

    plt.figure(figsize=(10, 8))
    plt.imshow(result_matrix, cmap='hot', interpolation='nearest', aspect='auto')
    plt.colorbar(label='Local Alignment Score (SW)')
    plt.title(f"Alignment Simulation: Influenza vs Covid\n(Windowing Method)")
    plt.xlabel("Covid-19 Windows")
    plt.ylabel("Influenza Windows")
    plt.show()