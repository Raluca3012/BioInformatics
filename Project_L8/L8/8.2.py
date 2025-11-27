#Download from NCBI a total of 3 bacterial genomes. Use these genomes as an 
#input for your app. Also, modify your app in order to be able to handle the
#amount of information from these genomes.
#Your app must detect transposable elements and the results from the output 
#must showw their position and length.
#In this case, the inverted repeats are unknown, ulnike the previous assignment, 
#in which they were known.
#Note! The following cases must be taken into consideration:
    #1. Transposons involving
    #2. Transposons overlapping
#NOTE2! The min size of the inverted digits must be of 4bases, and the max.
#of the inverted digits must be 6bases 

import matplotlib.pyplot as plt
import matplotlib.patches as patches

MIN_IR = 4
MAX_IR = 6
MIN_BODY = 20
MAX_BODY = 200
DR_LEN = 3

def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def read_fasta_sequence(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

def find_transposons_in_genome(seq):
    n = len(seq)
    results = []

    for L in range(MIN_IR, MAX_IR + 1):
        pos = {}
        for i in range(n - L + 1):
            k = seq[i:i+L]
            pos.setdefault(k, []).append(i)

        for left, left_list in pos.items():
            rc = revcomp(left)
            if rc not in pos:
                continue
            right_list = pos[rc]

            for i in left_list:
                min_j = i + L + MIN_BODY
                max_j = i + L + MAX_BODY
                if min_j >= n:
                    continue

                for j in right_list:
                    if j <= min_j:
                        continue
                    if j > max_j:
                        break
                    if j >= n - L - DR_LEN:
                        break
                    if i - DR_LEN < 0:
                        continue

                    DR_left = seq[i-DR_LEN:i]
                    DR_right = seq[j+L:j+L+DR_LEN]

                    if DR_left == DR_right:
                        start = i - DR_LEN
                        end = j + L + DR_LEN
                        results.append((start, end, L))

    results.sort()
    final = []
    seen = set()

    for s, e, L in results:
        if (s, e) not in seen:
            seen.add((s, e))
            final.append((s, e, L))

    return final

def find_overlaps(tes):
    pairs = []
    for i in range(len(tes)):
        s1, e1, _ = tes[i]
        for j in range(i+1, len(tes)):
            s2, e2, _ = tes[j]
            if s2 < e1:
                pairs.append((i, j))
    return pairs

def nucleotide_distribution(seq, title):
    counts = {b: seq.count(b) for b in "ACGT"}
    plt.figure(figsize=(6,4))
    plt.bar(counts.keys(), counts.values())
    plt.title("Nucleotide Distribution " + title)
    plt.xlabel("Base")
    plt.ylabel("Count")
    plt.grid(axis='y')
    plt.show()

def transposon_lengths(tes, title):
    if not tes:
        return
    labs = [f"TE{i+1}" for i in range(len(tes))]
    lens = [e - s for s, e, _ in tes]
    plt.figure(figsize=(7,4))
    plt.bar(labs, lens)
    plt.title("Transposon Lengths " + title)
    plt.xlabel("Transposon")
    plt.ylabel("Length (bp)")
    plt.grid(axis='y')
    plt.show()

def plot_transposon_structure(DR, IR_left, body, IR_right, title):
    fig, ax = plt.subplots(figsize=(12,3))
    x = 0
    DL = len(DR)
    IL = len(IR_left)
    BL = len(body)
    total = DL + IL + BL + IL + DL

    ax.text(x, 0.6, DR, color="red", fontsize=12)
    x += DL

    ax.add_patch(patches.Rectangle((x, 0.25), IL, 0.5, hatch='///', fill=False))
    ax.text(x, 0.05, "IR left")
    x += IL

    ax.add_patch(patches.Rectangle((x, 0.25), BL, 0.5, color="#b57edc"))
    ax.text(x + BL/2, 0.55, "Body", ha="center")
    x += BL

    ax.add_patch(patches.Rectangle((x, 0.25), IL, 0.5, hatch='///', fill=False))
    ax.text(x, 0.05, "IR right")
    x += IL

    ax.text(x, 0.6, DR, color="red", fontsize=12)

    ax.set_xlim(0, total + 2)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.title("Transposon Structure " + title)
    plt.show()

def summarize_genome(path):
    seq = read_fasta_sequence(path)
    print("\n=== Genome:", path, "===")
    print("Length:", len(seq))

    tes = find_transposons_in_genome(seq)
    print("Detected transposons:", len(tes))

    for i, (s, e, L) in enumerate(tes):
        print(f"TE{i+1}: start={s}, end={e}, length={e-s}, IR_len={L}")

    overlaps = find_overlaps(tes)
    print("Overlaps:", overlaps)

    nucleotide_distribution(seq, "for " + path)
    transposon_lengths(tes, "for " + path)

    if tes:
        s, e, L = tes[0]
        DR = seq[s:s+DR_LEN]
        IR_left = seq[s+DR_LEN:s+DR_LEN+L]
        IR_right = seq[e-DR_LEN-L:e-DR_LEN]
        body = seq[s+DR_LEN+L:e-DR_LEN-L]
        plot_transposon_structure(DR, IR_left, body, IR_right, "from " + path)

genomes = [
    "genome1.fasta",
    "genome2.fasta",
    "genome3.fasta"
]

for g in genomes:
    summarize_genome(g)

