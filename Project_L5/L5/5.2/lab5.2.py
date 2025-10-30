#Search and download from NCBI or other websites 10 viral genomes
#Use each of these genomes in order to take samples and make a set of samples
#for each viruse
#A) Measure the time (ms) of assembly for each set
#B) Measure their overall C+G percentage
#C) Plot a chart in which the coordinates of the points are as follows
# on the Y-axis, the point will take the time (ms) and on the X-axis the point
# will take the overall C+G percentage
#D) Make a text file in which you explain the differences between the positions 
# of the points.

import time
import random
import matplotlib.pyplot as plt

def read_fasta(path):
    seq = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()

def get_random_samples(S, num_samples=2000, min_len=100, max_len=150, seed=42):
    random.seed(seed)
    samples = []
    for _ in range(num_samples):
        start = random.randint(0, len(S) - max_len)
        frag_len = random.randint(min_len, max_len)
        frag = S[start:start + frag_len]
        samples.append(frag)
    return samples

def overlap_suffix_prefix(a, b, min_overlap=10):
    max_olen = min(len(a), len(b))
    for olen in range(max_olen, min_overlap - 1, -1):
        if a[-olen:] == b[:olen]:
            return olen, a + b[olen:]
    return 0, None

def greedy_assemble(reads, min_overlap=10):
    if not reads:
        return ""
    contig = reads.pop(0)
    changed = True
    while changed:
        changed = False
        for i, r in enumerate(reads):
            olen, merged = overlap_suffix_prefix(contig, r, min_overlap)
            if olen > 0:
                contig = merged
                reads.pop(i)
                changed = True
                break
            olen, merged = overlap_suffix_prefix(r, contig, min_overlap)
            if olen > 0:
                contig = merged
                reads.pop(i)
                changed = True
                break
    return contig

def cg_percentage(S):
    if not S:
        return 0
    cg = sum(1 for c in S if c in ('C', 'G'))
    return cg / len(S) * 100

if __name__ == "__main__":
    genomes = {
        "SARS-CoV-2": "covid.fasta",
        "Influenza A": "influentzaA.fasta",
        "Zika": "zika.fasta",
        "Dengue": "dengue.fasta",
        "Ebola": "ebola.fasta",
        "West Nile": "westNile.fasta",
        "MERS-CoV": "MERSCoV.fasta",
        "Papilloma": "papilloma.fasta",
        "Rabies": "rabies.fasta",
        "Norwalk": "norWalk.fasta"
    }

    times = []
    cg_perc = []
    labels = []

    for name, file in genomes.items():
        S = read_fasta(file)
        samples = get_random_samples(S, num_samples=2000, min_len=100, max_len=150, seed=7)
        start = time.time()
        contig = greedy_assemble(list(samples), min_overlap=10)
        elapsed = (time.time() - start) * 1000
        cg = cg_percentage(S)
        times.append(elapsed)
        cg_perc.append(cg)
        labels.append(name)
        print(f"{name}: {elapsed:.2f} ms, CG%={cg:.2f}")

    plt.figure(figsize=(10, 6))
    colors = plt.cm.tab10.colors

    for i, label in enumerate(labels):
        plt.scatter(cg_perc[i], times[i], color=colors[i % len(colors)], label=label, s=80)

    plt.xlabel("C+G percentage (%)")
    plt.ylabel("Assembly time (ms)")
    plt.title("Assembly Time vs C+G% for Viral Genomes")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
