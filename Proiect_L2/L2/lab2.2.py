#Find in sequence S only the dinucleotides and trinucleotides that exists, 
#without the use of brute force. In order to achieve the results one must verify 
#this combinations starting from the beg. of the seq.
def percentage(S, length):
    freq = {}
    total = len(S) - length + 1

    for i in range(total):
        part = S[i:i+length]
        if part in freq:
            freq[part] += 1
        else:
            freq[part] = 1

    for part in freq:
        freq[part] = round((freq[part] / total) * 100, 2)

    return freq


S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

percentage_din = percentage(S, 2)
print("Percentages for dinucleotides:")
for i in (percentage_din):
    print(f"{i}: {percentage_din[i]}%")

percentage_trin = percentage(S, 3)
print("Percentages for trinucleotides:")
for i in (percentage_trin):
    print(f"{i}: {percentage_trin[i]}%")
