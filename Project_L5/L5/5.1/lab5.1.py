#Take an arbitrary DNA sequence from NCBI, between 1000 and 3000 nucleotides
#A) Take 2000 random samples from this of about 100-150 bases
#B) Store these samples in an array variable/list
#C) Rebuild the original DNA seq by using these random samples
#D) What will be the main problem with your algorithm? Describe;

#Note what kind of structure inside the original seq may create different computation
#Note the samples must be aligned starting with the min of 10 positions in order to
#avoid random matching

import random
import matplotlib.pyplot as plt

def read_fasta(path):
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq.append(line)
    return "".join(seq).upper()

def get_random_samples(S, num_samples=2000, min_len=100, max_len=150, seed=42):
    random.seed(seed)
    samples = []
    if len(S) < max_len:
        raise ValueError("The seq. is too short to extract 1000-150 bases.")
    for _ in range(num_samples):
        start = random.randint(0, len(S) - max_len)
        frag_len = random.randint(min_len, max_len)
        frag = S[start:start + frag_len]
        samples.append((start, frag))
    return samples

def overlap_suffix_prefix(a, b, min_overlap=10):
    max_olen = min(len(a), len(b))
    for olen in range(max_olen, min_overlap - 1, -1):
        if a[-olen:] == b[:olen]:
            return olen, a + b[olen:]
    return 0, None

def greedy_assemble(reads, min_overlap=10, max_passes=2000):
    if not reads:
        return "", 0
    contig = reads.pop(0)
    merges = 0
    for _ in range(max_passes):
        changed = False
        i = 0
        while i < len(reads):
            olen, merged = overlap_suffix_prefix(contig, reads[i], min_overlap)
            if olen > 0:
                contig = merged
                reads.pop(i)
                merges += 1
                changed = True
                i = 0
            else:
                i += 1
        i = 0
        while i < len(reads):
            olen, merged = overlap_suffix_prefix(reads[i], contig, min_overlap)
            if olen > 0:
                contig = merged
                reads.pop(i)
                merges += 1
                changed = True
                i = 0
            else:
                i += 1
        if not changed:
            break
    return contig, merges

def reconstruct_from_samples(samples, original_len, min_overlap=10):
    reads = [frag for (_, frag) in samples]
    if not reads:
        return "N" * original_len
    contig, merges = greedy_assemble(reads, min_overlap=min_overlap)
    if len(contig) < original_len:
        contig = contig + ("N" * (original_len - len(contig)))
    elif len(contig) > original_len:
        contig = contig[:original_len]
    return contig

def print_main_problem_notes():
    print("\n[D] Main problem with the algorithm:")
    print("- Repetitive sequences can cause false overlaps.")
    print("- Coverage may be uneven: some regions may not appear in the samples → 'N' in reconstruction.")
    print("- Greedy assembly does not guarantee the global optimum; it may stop at an incorrect arrangement.")
    print("- Even with a 10-base overlap threshold, random matches can still occur in repetitive regions.")

if __name__ == "__main__":
    fasta_path = "lab5.fasta"
    S = read_fasta(fasta_path)
    if not (1000 <= len(S) <= 3000):
        print("The seq. has not between 1000-3000 nucleotides.")
    samples = get_random_samples(S, num_samples=2000, min_len=100, max_len=150, seed=7)
    rebuilt = reconstruct_from_samples(samples, original_len=len(S), min_overlap=10)
    matches = sum(1 for a, b in zip(S, rebuilt) if a == b)
    identity = matches / len(S) * 100 if S else 0.0
    print(f"Identity with original: {identity:.2f}%")
    print("Original - first 120: ", S[:120])
    print("Rebuilt - first 120: ", rebuilt[:120])
    print_main_problem_notes()

    plt.figure(figsize=(12, 6))
    y_offset = 0

    plt.hlines(y=y_offset, xmin=0, xmax=len(S), color='black', linewidth=2, label='Original sequence')
    y_offset -= 1

    subset = random.sample(samples, 60)
    for i, (start, frag) in enumerate(subset):
        plt.hlines(y=y_offset - i, xmin=start, xmax=start + len(frag), color='blue', linewidth=3)
        plt.text(start, y_offset - i + 0.15, f"S{i+1}", fontsize=8, color='blue')

    plt.hlines(y=y_offset - len(subset) - 1, xmin=0, xmax=len(rebuilt), color='red', linewidth=2, label='Reconstructed sequence (C)')

    plt.title("Alignment Scheme of Random DNA Fragments")
    plt.xlabel("Position along the sequence")
    plt.ylabel("Fragments (S₁ ... Sₙ and C)")
    plt.legend()
    plt.xlim(0, len(S))
    plt.ylim(y_offset - len(subset) - 2, 2)
    plt.tight_layout()
    plt.show()
