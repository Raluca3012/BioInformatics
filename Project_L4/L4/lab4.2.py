# Download from NCBI the FASTA files containing the COVID-19 genome, and the 
# influenza genome. Use AI to compare the codon frequencies between the 2 genomes.
# A) Make a chart that shows the top 10 most frequent for COVID-19
# B) Make a chart that shows the top 10 most frequent for influenza
# C) Compare the 2 results and show the most freq. codors above 2 genoms
# D) Show in the output of the console the top 3 aminoacids for each genome
# E) Formulate a prompt for AI such that the 3 aminoacids are used to ask the AI
# which foods contain less those aminoacids

from collections import Counter
import matplotlib.pyplot as plt

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    seq = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return seq.upper()

def codon_frequencies(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]
    return Counter(codons)

def codon_to_aminoacid(codon):
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }
    return table.get(codon, '?')

def top_codons(freqs, n=10):
    return freqs.most_common(n)

def plot_top_codons(freqs, title):
    plt.figure()
    codons, counts = zip(*freqs)
    total = sum(counts)
    percentages = [round((c / total) * 100, 2) for c in counts]

    plt.bar(codons, percentages, color='royalblue', edgecolor='black')
    plt.title(title)
    plt.xlabel("Codons")
    plt.ylabel("Frequency (%)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    for i, v in enumerate(percentages):
        plt.text(i, v + max(percentages)*0.01, f"{v}%", ha='center', fontsize=8)

    plt.tight_layout()
    plt.show()

def aminoacid_frequencies(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    aminoacids = [codon_to_aminoacid(c) for c in codons if len(c) == 3]
    return Counter(aminoacids)

covid_seq = read_fasta("covid.fasta")
flu_seq = read_fasta("influenza.fna")

covid_freqs = codon_frequencies(covid_seq)
flu_freqs = codon_frequencies(flu_seq)

# A) Top 10 codons COVID
top10_covid = top_codons(covid_freqs)
print("Top 10 codons COVID:")
for c, v in top10_covid:
    print(f"{c}: {v}")
plot_top_codons(top10_covid, "Top 10 codons - COVID-19")

# B) Top 10 codons Influenza
top10_flu = top_codons(flu_freqs)
print("\nTop 10 codons Influenza:")
for c, v in top10_flu:
    print(f"{c}: {v}")
plot_top_codons(top10_flu, "Top 10 codons - Influenza")

# C) Compare the 2 results and show the most frequent codons across both genomes

common_codons = set([c for c, _ in top10_covid]) & set([c for c, _ in top10_flu])

print("\nCommon frequent codons between COVID and Influenza:")
print(", ".join(common_codons) if common_codons else "None")

# dacă există codoni comuni, îi comparăm grafic
if common_codons:
    covid_common = {c: covid_freqs[c] for c in common_codons}
    flu_common = {c: flu_freqs[c] for c in common_codons}

    labels = list(common_codons)
    covid_vals = [covid_common[c] for c in labels]
    flu_vals = [flu_common[c] for c in labels]

    total_covid = sum(covid_freqs.values())
    total_flu = sum(flu_freqs.values())
    covid_percent = [round(v / total_covid * 100, 2) for v in covid_vals]
    flu_percent = [round(v / total_flu * 100, 2) for v in flu_vals]

    x = range(len(labels))
    plt.figure()
    plt.bar(x, covid_percent, width=0.4, label="COVID-19", color='royalblue', edgecolor='black')
    plt.bar([i + 0.4 for i in x], flu_percent, width=0.4, label="Influenza", color='orange', edgecolor='black')
    plt.xticks([i + 0.2 for i in x], labels)
    plt.ylabel("Frequency (%)")
    plt.title("Common Frequent Codons - COVID-19 vs Influenza")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()


# D) Top 3 amino acids for each genome
covid_aas = aminoacid_frequencies(covid_seq)
flu_aas = aminoacid_frequencies(flu_seq)

top3_covid_aa = covid_aas.most_common(3)
top3_flu_aa = flu_aas.most_common(3)

print("\nTop 3 amino acids (COVID):", top3_covid_aa)
print("Top 3 amino acids (Influenza):", top3_flu_aa)


