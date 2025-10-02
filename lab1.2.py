def relative_frequency(S):
    freq = {}
    n = len(S)
    for char in S:
        freq[char] = freq.get(char, 0) + 1
    for char in freq:
        freq[char]/=n
    return freq


S = input('Input a sequence of 20 letters: ')
print('The relative freqs. are: ', relative_frequency(S))