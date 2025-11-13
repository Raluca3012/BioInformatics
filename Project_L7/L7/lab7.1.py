def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()


def find_repeats(seq, min_len, max_len, min_reps):
    repeats = {}
    for L in range(min_len, max_len + 1):
        seen = {}
        for i in range(len(seq) - L + 1):
            frag = seq[i:i + L]
            seen[frag] = seen.get(frag, 0) + 1
        for frag, count in seen.items():
            if count >= min_reps:
                repeats[(frag, L)] = count
    return repeats


def print_repeats(reps):
    with open("repeats_output.txt", "w") as f:
        if not reps:
            f.write("No repetitive seq. found.\n")
            return

        for (frag, L), cnt in sorted(reps.items(), key=lambda x: (-x[1], x[0][1])):
            line = f"{frag} (len={L})  ->  {cnt} times\n"
            f.write(line)

seq = read_fasta("covid.fasta")[:3000]

if not (1000 <= len(seq) <= 3000):
    print("WARNING: Seq. length is not within 1000–3000 bp.")

reps = find_repeats(seq, 3, 6, 2)
print_repeats(reps)
