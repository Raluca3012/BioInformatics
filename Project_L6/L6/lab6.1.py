import random
import math
import matplotlib.pyplot as plt

def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()

def random_fragments(seq, n=10, min_len=100, max_len=3000):
    fragments = []
    L = len(seq)
    for _ in range(n):
        size = random.randint(min_len, min(max_len, L))
        start = random.randint(0, L - size)
        fragment = seq[start:start + size]
        fragments.append(fragment)
    return fragments

def ascii_gel(fragments):
    sizes = sorted([len(f) for f in fragments], reverse=True)
    max_bp = max(sizes)
    print("ASCII GEL REPRESENTATION\n")
    for bp in sizes:
        pos = int((1 - bp / max_bp) * 40)
        print(f"{bp:4} bp | " + " " * pos + "█")

def plot_gel(fragments):
    sizes = sorted([len(f) for f in fragments], reverse=True)
    max_bp = max(sizes)
    min_bp = min(sizes)

    def gel_pos(bp):
        return (math.log10(max_bp) - math.log10(bp)) / (math.log10(max_bp) - math.log10(min_bp))

    fig, ax = plt.subplots(figsize=(4,7), dpi=150)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


    for bp in sizes:
        y = 0.9 - gel_pos(bp) * 0.8
        ax.hlines(y, 0.32, 0.68, linewidth=3)

    for mark in [3000, 2500, 1000, 500, 150]:
        if min_bp <= mark <= max_bp:
            y = 0.9 - gel_pos(mark) * 0.8
            ax.text(0.1, y, f"{mark} bp -", color="black", fontsize=9)

    ax.set_title("Simulated Gel Electrophoresis", pad=8)
    plt.tight_layout()
    plt.show()


full_seq = read_fasta("ebola.fasta")
print("Full sequence length:", len(full_seq), "bp")

seq = full_seq[:3000]
print("Selected sequence length:", len(seq), "bp")

fragments = random_fragments(seq)
print("Fragments stored in array (bp):")
print([len(f) for f in fragments])

ascii_gel(fragments)
plot_gel(fragments)
