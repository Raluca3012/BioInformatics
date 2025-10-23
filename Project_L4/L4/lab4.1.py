#Implement an application that converts the coding region of a gene into an
#aminoacid sequence. Use the genetic code table from Moodle.
#Take a seq of DNA and convert it into ARN seq:
    # SDNA = ...ATG......TAA...
    # SRNA = ...AUG......UAA... 
    
def dna_to_rna(S):
    return S.replace('T', 'U')

def translate_rna(S):
    codon_table = {
        'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
        'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
        'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
        'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
        'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
        'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'UAU':'Y', 'UAC':'Y', 'UAA':'*', 'UAG':'*',
        'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
        'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
        'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'UGU':'C', 'UGC':'C', 'UGA':'*', 'UGG':'W',
        'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
        'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }

    start = S.find('AUG')
    if start == -1:
        return "No START codon found"
    
    protein = ''
    for i in range(start, len(S) - 2, 3):
        codon = S[i:i+3]
        aa = codon_table.get(codon, '')
        if aa == '*':
            break
        protein += aa
    return protein

S_DNA = "AATGGCTTCGACG"
S_RNA = dna_to_rna(S_DNA)
S_PROT = translate_rna(S_RNA)

print("DNA sequence: ", S_DNA)
print("RNA sequence: ", S_RNA)
print("Aminoacid sequence: ", S_PROT)
