import os
import matplotlib.pyplot as plt

W = 30

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq.append(line.upper())
    return "".join(seq)

def CG_content(win):
    c = win.count("C")
    g = win.count("G")
    return round((c+g)/len(win)*100, 2)

def IC_promkappa(win):
    A = win
    L = len(A)
    shifts = L - 1
    total = 0.0

    for u in range(1, shifts+1):
        seg = L - u
        matches = 0
        for i in range(seg):
            if A[i] == A[i+u]:
                matches += 1
        total += (matches/seg)*100.0

    return round(total/shifts, 2)

def sliding(seq, w):
    return [seq[i:i+w] for i in range(len(seq)-w+1)]

def compute_pattern(seq, w):
    wins = sliding(seq, w)
    xs = [CG_content(wd) for wd in wins]
    ys = [IC_promkappa(wd) for wd in wins]
    cx = sum(xs)/len(xs)
    cy = sum(ys)/len(ys)
    return xs, ys, cx, cy

def main():

    root = os.getcwd()
    folders = ["influenza", "covid"]    
    centers = []
    plt.figure(figsize=(14,6))
    plt.title("ODS Pattern for Influenza & COVID Genomes")
    plt.xlabel("Mean C+G%")
    plt.ylabel("Kappa IC")
    plt.grid(True)
    for folder in folders:
        folder_path = os.path.join(root, folder)
        for file in os.listdir(folder_path):
            if not file.lower().endswith(".fasta"):
                continue 
            fpath = os.path.join(folder_path, file)
            seq = read_fasta(fpath)
            xs, ys, cx, cy = compute_pattern(seq, W)
            centers.append((cx, cy, file.replace(".fasta","")))
            plt.scatter(xs, ys, alpha=0.45, s=8, label=file)

    plt.legend(fontsize=6)
    plt.tight_layout()
    plt.show()
    plt.figure(figsize=(10,6))
    plt.title("Centers of Weight for Each Genome")
    plt.xlabel("Avg C+G%")
    plt.ylabel("Avg IC")
    plt.grid(True)

    for cx, cy, name in centers:
        plt.scatter(cx, cy, s=90)
        plt.text(cx+0.1, cy+0.1, name, fontsize=7)

    plt.show()


if __name__ == "__main__":
    main()
