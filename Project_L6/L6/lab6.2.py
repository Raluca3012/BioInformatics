# Download 10 influenza genomes and apply the electrophoresis gel selection 
#on each of them. Make a comparison between the electrophoresis gel simulations
# and show which of the influenza genome show the most DNA seq.
# You can plot them in the same graph, but also separatelly because the lines may overlap
# As the main restriction enzimes please use ECOR1

# Electrophoresis gel simulation for 10 influenza genomes with EcoRI
# Reads 10 FASTA files, digests with EcoRI, keeps fragments in [100, 3000] bp,
# plots a combined gel (all lanes) and individual gels,
# reports which genome shows the most DNA (sum of fragment lengths in window) and most bands.

import os
import math
import matplotlib.pyplot as plt

def read_fasta(filename):
    seq = []
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()

def ecoRI_digest(seq):
    site = "GAATTC"
    parts = []
    i = 0
    while True:
        j = seq.find(site, i)
        if j == -1:
            parts.append(seq[i:])
            break
        parts.append(seq[i:j+1])
        i = j+1
    return parts

def filter_window(fragments, min_bp=100, max_bp=3000):
    return [f for f in fragments if min_bp <= len(f) <= max_bp]

def gel_pos_mapper(sizes, floor_min=100, ceil_max=3000):
    if not sizes:
        return lambda bp: 0.5
    min_bp = max(floor_min, min(sizes))
    max_bp = max(ceil_max, max(sizes))
    lmin = math.log10(min_bp)
    lmax = math.log10(max_bp)
    def pos(bp):
        x = max(min_bp, min(bp, max_bp))
        return (math.log10(max_bp) - math.log10(x)) / (lmax - lmin + 1e-12)
    return pos

def plot_gel_combined(name_to_sizes, markers=(3000,1000,500)):
    labels = list(name_to_sizes.keys())
    lanes = [name_to_sizes[k] for k in labels]
    all_sizes = [s for sizes in lanes for s in sizes] or [100,3000]
    gpos = gel_pos_mapper(all_sizes)

    COLORS = [
        "red", "blue", "green", "purple", "orange",
        "brown", "cyan", "magenta", "gold", "black"
    ]

    fig, ax = plt.subplots(figsize=(11,7), dpi=130)
    ax.set_xlim(0, len(labels)+1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    for idx, sizes in enumerate(lanes, start=1):
        color = COLORS[(idx-1) % len(COLORS)]
        ax.add_patch(plt.Rectangle((idx-0.15, 0.95), 0.30, 0.03, fill=False, edgecolor=color, linewidth=2))
        for bp in sorted(sizes, reverse=True):
            y = 0.92 - gpos(bp) * 0.80
            ax.hlines(y, idx-0.12, idx+0.12, color=color, linewidth=3)
        ax.text(idx, 0.02, labels[idx-1], ha="center", va="bottom", fontsize=9, rotation=45, color=color)

    for m in markers:
        y = 0.92 - gpos(m) * 0.80
        ax.hlines(y, 0.4, 0.8, linewidth=3, color="black")
        ax.text(0.25, y, f"{m} bp", va="center", fontsize=9)

    ax.set_title("EcoRI", pad=8)
    plt.tight_layout()
    plt.show()


files = [
    ("ebola.fasta", "Ebola"),
    ("covid.fasta", "COVID-19"),
    ("dengue.fasta", "Dengue"),
    ("influentzaA.fasta", "Influenza A"),
    ("papilloma.fasta", "Papillomavirus"),
    ("rabies.fasta", "Rabies"),
    ("westNile.fasta", "West Nile virus"),
    ("zika.fasta", "Zika"),
    ("MERSCoV.fasta", "MERS-CoV"),
    ("norWalk.fasta", "Norwalk"),
]

existing = [(p, n) for (p, n) in files if os.path.exists(p)]
if len(existing) < 1:
    raise FileNotFoundError("Fișierele .fasta nu au fost găsite în folderul curent.")

name_to_sizes = []
report = []

for path, nice_name in existing:
    seq = read_fasta(path)
    frags = ecoRI_digest(seq)
    frags_sel = filter_window(frags, 100, 3000)
    sizes = [len(f) for f in frags_sel]
    name_to_sizes.append((nice_name, sizes))
    total_bp = sum(sizes)
    report.append((nice_name, total_bp, len(sizes), max(sizes) if sizes else 0))

name_to_sizes_dict = dict(name_to_sizes)
plot_gel_combined(name_to_sizes_dict)

report.sort(key=lambda x: (-x[1], -x[2], -x[3]))

for name, total_bp, band_count, max_frag in report:
    print(f"{name}: total_bp={total_bp}, bands={band_count}, max_frag={max_frag}")

if report:
    print("\nMost DNA on gel:", report[0][0])
