#https://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/Bacteroides_4_3_47FAA_uid32443/NZ_ACDR00000000.scaffold.faa.tgz

import tkinter as tk
from tkinter import ttk, messagebox
import tarfile
import os

def detect_alphabet(S):
    return set(S)

def relative_frequencies(S):
    freq = {}
    total = len(S)
    for char in S:
        freq[char] = freq.get(char, 0) + 1
    for char in freq:
        freq[char] /= total
    return freq

def read_fasta(file_path):
    #Read a FASTA file and return a list of (header, sequence)
    sequences = []
    with open(file_path, "r") as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # new sequence
                if header and seq_lines:
                    sequences.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:  # last sequence
            sequences.append((header, "".join(seq_lines)))
    return sequences

def analyze_fasta():
    try:
        tgz_path = r"C:\Users\rdoro\OneDrive\Desktop\Anul_IV\Bioinformatics\lab1.3\NZ_ACDR00000000.scaffold.faa.tgz"
        extract_dir = os.path.dirname(tgz_path)
        
        progress_bar["value"] = 10
        root.update_idletasks()
        
        with tarfile.open(tgz_path, "r:gz") as tar:
            tar.extractall(path=extract_dir)

        progress_bar["value"] = 40
        root.update_idletasks()

        faa_files = [f for f in os.listdir(extract_dir) if f.endswith(".faa")]
        if not faa_files:
            messagebox.showerror("Error", "No .faa files found after extraction!")
            return
        faa_path = os.path.join(extract_dir, faa_files[0])

        sequences = read_fasta(faa_path)
        if not sequences:
            messagebox.showerror("Error", "No sequences found in the file!")
            return

        progress_bar["value"] = 70
        root.update_idletasks()

        result = f"Found {len(sequences)} sequences in the file.\n\n"
        for i, (header, seq) in enumerate(sequences[:3], start=1):  
            alphabet = detect_alphabet(seq)
            freqs = relative_frequencies(seq)
            result += f"Sequence {i}:\n{header}\nLength: {len(seq)}\n"
            result += f"Alphabet: {alphabet}\nRelative frequencies:\n"
            for k, v in sorted(freqs.items()):
                result += f"  {k}: {v:.3f}\n"
            result += "\n"

        global_seq = "".join(seq for _, seq in sequences)
        alphabet = detect_alphabet(global_seq)
        freqs = relative_frequencies(global_seq)
        result += "\n=== Global Summary for All Sequences ===\n"
        result += f"Total alphabet: {alphabet}\nGlobal relative frequencies:\n"
        for k, v in sorted(freqs.items()):
            result += f"  {k}: {v:.3f}\n"

        progress_bar["value"] = 100
        root.update_idletasks()

        messagebox.showinfo("FASTA Analysis Result", result)

    except Exception as e:
        messagebox.showerror("Error", str(e))

root = tk.Tk()
root.title("FASTA Analyzer (.tgz direct)")

frame = tk.Frame(root)
frame.pack(pady=20, padx=20)

button = tk.Button(frame, text="Analyze FASTA File", command=analyze_fasta, font=("Arial", 14))
button.pack(pady=10)

progress_bar = ttk.Progressbar(frame, orient="horizontal", length=300, mode="determinate")
progress_bar.pack(pady=10)

root.mainloop()
