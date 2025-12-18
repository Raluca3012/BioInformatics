import math

BASES = "ACGT"
NULL_P = 0.25

MOTIFS = ["GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
          "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT",]

def print_motif_table(motifs):
    L = len(motifs[0])
    print("    " + "  ".join(str(i+1) for i in range(L)))
    print("   " + "---"*L)
    for i, m in enumerate(motifs, 1):
        print(f"{i:2d} | " + "  ".join(m))
    print()

def build_count_matrix(motifs):
    L = len(motifs[0])
    return {
        b: [sum(m[i] == b for m in motifs) for i in range(L)]
        for b in BASES
    }

def build_relative_freq_matrix(counts, n):
    return {
        b: [counts[b][i] / n for i in range(len(counts[b]))]
        for b in BASES
    }

def build_log_likelihood_matrix(freqs):
    return {
        b: [
            0.0 if freqs[b][i] == 0 else math.log(freqs[b][i] / NULL_P)
            for i in range(len(freqs[b]))
        ]
        for b in BASES
    }

def print_matrix(title, matrix, fmt="{:6.2f}"):
    L = len(next(iter(matrix.values())))
    print(title)
    print("    " + " ".join(f"{i+1:>6}" for i in range(L)))
    for b in BASES:
        row = " ".join(fmt.format(matrix[b][i]) for i in range(L))
        print(f"{b} : {row}")
    print()

print("Known motif sequences (aligned):")
print_motif_table(MOTIFS)

counts = build_count_matrix(MOTIFS)
freqs = build_relative_freq_matrix(counts, len(MOTIFS))
logll = build_log_likelihood_matrix(freqs)

print_matrix("Count matrix", counts, fmt="{:6d}")
print_matrix("Relative frequencies matrix", freqs, fmt="{:6.2f}")
print_matrix("Log-likelihood matrix", logll, fmt="{:6.2f}")


S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

def score_window(window, logll):
    score = 0.0
    for i, base in enumerate(window):
        v = logll[base][i]
        if v == 0.0:
            return 0.0
        score += v
    return score

L = len(MOTIFS[0])
threshold = min(score_window(m, logll) for m in MOTIFS)

print("Threshold:", round(threshold, 3))
print("Index  Window        Score")

for i in range(len(S) - L + 1):
    w = S[i:i+L]
    sc = score_window(w, logll)
    mark = " <-- signal" if sc >= threshold else ""
    print(f"{i:>5}  {w}  {sc:>8.3f}{mark}")
