import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def make_transposon():
    DR = random_dna(3)
    IR_left = random_dna(4)
    IR_right = revcomp(IR_left)
    body = random_dna(random.randint(15, 25))
    TE = DR + IR_left + body + IR_right + DR
    return TE, DR, IR_left

def insert_TE(host, TE, pos):
    return host[:pos] + TE + host[pos:]

def generate_sequence():
    host = random_dna(random.randint(200, 350))
    TEs = []
    positions = []
    for _ in range(3):
        TE, DR, IR_left = make_transposon()
        TEs.append((TE, DR, IR_left))
    pos1 = random.randint(30, 80)
    host = insert_TE(host, TEs[0][0], pos1)
    positions.append(pos1)
    pos2 = random.randint(120, 180)
    host = insert_TE(host, TEs[1][0], pos2)
    positions.append(pos2)
    pos3 = pos2 + random.randint(5, 15)
    host = insert_TE(host, TEs[2][0], pos3)
    positions.append(pos3)
    return host, TEs, positions

def detect_transposons(sequence, TEs):
    results = []
    for TE, DR, IR_left in TEs:
        IR_right = revcomp(IR_left)
        len_IR = len(IR_left)
        len_DR = len(DR)
        for i in range(len(sequence)):
            if sequence[i:i+len_DR] == DR and sequence[i+len_DR:i+len_DR+len_IR] == IR_left:
                for j in range(i+len_DR+len_IR, len(sequence)-len_IR):
                    if sequence[j:j+len_IR] == IR_right:
                        rightDR = sequence[j+len_IR : j+len_IR+len_DR]
                        if rightDR == DR:
                            results.append((i, j+len_IR+len_DR))
    return results

def plot_transposon_structure(DR, IR_left, body, IR_right):
    fig, ax = plt.subplots(figsize=(12, 3))
    x = 0
    DR_len = len(DR)
    IR_len = len(IR_left)
    body_len = len(body)
    total = DR_len + IR_len + body_len + IR_len + DR_len
    ax.text(x, 0.6, DR, color="red", fontsize=12)
    x += DR_len
    ax.add_patch(patches.Rectangle((x, 0.25), IR_len, 0.5, hatch='///', fill=False))
    ax.text(x, 0.1, "IR left", fontsize=10)
    x += IR_len
    ax.add_patch(patches.Rectangle((x, 0.25), body_len, 0.5, color="#b57edc"))
    ax.text(x + body_len/2, 0.55, "Transposon body", ha="center", fontsize=10)
    x += body_len
    ax.add_patch(patches.Rectangle((x, 0.25), IR_len, 0.5, hatch='///', fill=False))
    ax.text(x, 0.1, "IR right", fontsize=10)
    x += IR_len
    ax.text(x, 0.6, DR, color="red", fontsize=12)
    ax.set_xlim(-1, total + 2)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.title("Transposon Structure")
    plt.show()

host, TEs, real_positions = generate_sequence()
print("Artificial DNA sequence")
print(host)
print("\nReal TE insert positions:", real_positions)
found = detect_transposons(host, TEs)
print("\nDETECTED TRANSPOSONS")
lengths = []
for start, end in found:
    L = end - start
    lengths.append(L)
    print(f"Transposon found from {start} to {end} (length {L} bp)")
counts = {"A": 0, "C": 0, "G": 0, "T": 0}
for base in host:
    counts[base] += 1
plt.figure(figsize=(8, 5))
plt.bar(counts.keys(), counts.values(), color=["red", "blue", "green", "orange"])
plt.title("Nucleotide Distribution in Artificial DNA")
plt.xlabel("Nucleotide")
plt.ylabel("Frequency")
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.show()
if lengths:
    plt.figure(figsize=(8, 5))
    plt.bar([f"TE{i+1}" for i in range(len(lengths))], lengths, color="purple")
    plt.title("Lengths of Detected Transposable Elements")
    plt.xlabel("Transposon")
    plt.ylabel("Length (bases)")
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.show()
TE, DR, IR_left = TEs[0]
IR_right = revcomp(IR_left)
body = TE[len(DR) + len(IR_left) : -(len(DR) + len(IR_left))]
plot_transposon_structure(DR, IR_left, body, IR_right)
