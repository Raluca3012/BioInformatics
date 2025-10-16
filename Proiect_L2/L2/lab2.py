from itertools import product

def percentage(S, length):
    freq = {}
    total = len(S) - length + 1 
    for i in range (total):
        part = S[i:i+length]
        if part in freq:
            freq[part]+=1
        else:
            freq[part] = 1
            
    return freq

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
nucleotides = ['A', 'C', 'G', 'T']
dinucleotides = [''.join(p) for p in product(nucleotides, repeat = 2)]
trinucleotides = [''.join(p) for p in product(nucleotides, repeat = 3)]

print("Dinucleotides: ", dinucleotides)
print("Trinucleotides: ", trinucleotides)

percentage_din = percentage(S, 2)
print("Percentages for dinucleotides: ")
for i in sorted(percentage_din):
    print(f"{i}:  {percentage_din[i]}%")
    
percentage_trin = percentage(S, 3)
print("Percentages for trinucleotides: ")
for i in sorted(percentage_trin):
    print(f"{i}:  {percentage_trin[i]}%")