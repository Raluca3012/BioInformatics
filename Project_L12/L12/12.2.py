#Download 10 influenza genomes , adapt  your application from previous 
#assignment in order to scan each genome for possible motives . For each genome
#make a chart that shows the signal with most likely locations of real 
#functional motives 

import math
import os
import matplotlib.pyplot as plt

BASES = "ACGT"
NULL_P = 0.25
L = 9
plt.close('all')
MOTIFS = ["GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
          "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT",]

FASTA_DIR = "influenza"

def build_count_matrix(motifs):
    return {
        b: [sum(m[i] == b for m in motifs) for i in range(L)]
        for b in BASES
    }

def build_relative_freq_matrix(counts, n):
    return {
        b: [counts[b][i] / n for i in range(L)]
        for b in BASES
    }

def build_log_likelihood_matrix(freqs):
    return {
        b: [
            0.0 if freqs[b][i] == 0 else math.log(freqs[b][i] / NULL_P)
            for i in range(L)
        ]
        for b in BASES
    }

def score_window(window, logll):
    s = 0.0
    for i, base in enumerate(window):
        v = logll[base][i]
        if v == 0.0:
            return 0.0
        s += v
    return s

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

counts = build_count_matrix(MOTIFS)
freqs = build_relative_freq_matrix(counts, len(MOTIFS))
logll = build_log_likelihood_matrix(freqs)

threshold = min(score_window(m, logll) for m in MOTIFS)
print("Threshold:", round(threshold, 3))

for fname in sorted(os.listdir(FASTA_DIR)):
    if not fname.endswith(".fasta"):
        continue

    path = os.path.join(FASTA_DIR, fname)
    seq = read_fasta(path)

    scores = []
    for i in range(len(seq) - L + 1):
        w = seq[i:i+L]
        if any(c not in BASES for c in w):
            scores.append(0.0)
        else:
            scores.append(score_window(w, logll))

    best_idx = max(range(len(scores)), key=lambda i: scores[i])
    best_score = scores[best_idx]

    print(f"{fname}: best score = {best_score:.3f} at position {best_idx}")

    plt.figure()
    plt.plot(scores)   
    plt.axhline(threshold, color="red", linestyle="--", linewidth=1.5, label="Threshold")
    plt.xlabel("Sliding window start index")
    plt.ylabel("Log-likelihood score")
    plt.title(fname)
    plt.show()
