#Download 10 influenza genomes, for each genome plot on a chart the most
# frequent repetitions.

import matplotlib.pyplot as plt
import os

def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()

def find_repeats(seq, min_len=3, max_len=6, min_reps=2):
    repeats = {}
    for L in range(min_len, max_len + 1):
        seen = {}
        for i in range(len(seq) - L + 1):
            frag = seq[i:i + L]
            seen[frag] = seen.get(frag, 0) + 1
        for frag, c in seen.items():
            if c >= min_reps:
                repeats[(frag, L)] = c
    return repeats

files = [
    "A_RVA.fasta", "A_California.fasta", "A_NY.fasta", "A_Puerto_Rico.fasta", 
    "A_Rabies.fasta", "A_Victoria.fasta", "B_Brisbane.fasta",
    "C_AnnArbor.fasta", "A_Tennessee.fasta", "D_Oklahoma.fasta"
]

names = []
values = []

for f in files:
    if not os.path.exists(f):
        print("Missing:", f)
        continue

    seq = read_fasta(f)
    seq = seq[:3000]

    reps = find_repeats(seq, 3, 6, 2)

    if reps:
        best_pattern, best_count = max(reps.items(), key=lambda x: x[1])
        print(f, "->", best_pattern[0], "(len", best_pattern[1], ")", "appears", best_count, "times")
        names.append(f)
        values.append(best_count)
    else:
        print(f, "-> no repeats")
        names.append(f)
        values.append(0)

plt.figure(figsize=(12,6))
plt.bar(names, values)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Most frequent repetition count")
plt.title("Top repetition frequency in 10 influenza genomes")
plt.tight_layout()
plt.show()
