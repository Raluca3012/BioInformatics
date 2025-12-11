import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def needleman_wunsch(seq1, seq2, match_reward=1, mismatch_penalty=-1, gap_penalty=0):
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    for i in range(n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0][j] = j * gap_penalty

    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                match = score_matrix[i-1][j-1] + match_reward
            else:
                match = score_matrix[i-1][j-1] + mismatch_penalty
                
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            max_score = max(match, delete, insert)
            score_matrix[i][j] = max_score
            
            if max_score == match:
                traceback_matrix[i][j] = 1 
            elif max_score == delete:
                traceback_matrix[i][j] = 2 
            else:
                traceback_matrix[i][j] = 3
                
    align1 = ""
    align2 = ""
    i, j = n, m
    
    path_coords = [] 
    path_coords.append((j, i)) 
    
    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        is_match = (i > 0 and j > 0 and (seq1[i-1] == seq2[j-1]))
        diag_score = score_matrix[i-1][j-1] if (i > 0 and j > 0) else -float('inf')
        up_score = score_matrix[i-1][j] if i > 0 else -float('inf')
        left_score = score_matrix[i][j-1] if j > 0 else -float('inf')
        match_val = match_reward if is_match else mismatch_penalty

        if i > 0 and j > 0 and current_score == diag_score + match_val:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and current_score == up_score + gap_penalty:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1
            
        path_coords.append((j, i))

    align1 = align1[::-1]
    align2 = align2[::-1]
    
    return align1, align2, score_matrix, path_coords

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

GAP_PENALTY = 0     
MATCH_SCORE = 1   
MISMATCH_SCORE = -1  

aln1, aln2, mat, path = needleman_wunsch(S1, S2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY)

matches = sum(1 for a, b in zip(aln1, aln2) if a == b)
length = len(aln1)
similarity = (matches / length) * 100

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

colors = ["black", "red"]
cmap_heatmap = mcolors.LinearSegmentedColormap.from_list("black_red", colors)

im = axes[0].imshow(mat, cmap=cmap_heatmap, aspect='auto')
axes[0].set_title("Graphic representation of the alignment matrix")
axes[0].set_xlabel("Sequence 2 (S2)")
axes[0].set_ylabel("Sequence 1 (S1)")
plt.colorbar(im, ax=axes[0])

axes[1].set_title("Traceback path")
axes[1].set_xlim(0, len(S2))
axes[1].set_ylim(0, len(S1))
axes[1].set_xticks(np.arange(0, len(S2) + 1, 1))
axes[1].set_yticks(np.arange(0, len(S1) + 1, 1))
axes[1].grid(color='black', linestyle='-', linewidth=0.5)
axes[1].set_facecolor('#FFFFE0') 

for x, y in path:
    rect = plt.Rectangle((x, len(S1) - y - 1), 1, 1, facecolor='red', edgecolor='black')
    axes[1].add_patch(rect)

axes[1].set_xticklabels([''] + list(S2))
axes[1].set_yticklabels([''] + list(S1)[::-1])
axes[1].set_xlabel("S2")
axes[1].set_ylabel("S1")

plt.tight_layout()
plt.show()

print(f"Alignment Results:")
print(f"S1: {aln1}")
print(f"    " + "".join(['|' if a==b else ' ' for a,b in zip(aln1, aln2)]))
print(f"S2: {aln2}")
print("-" * 30)
print(f"Matches    = {matches}")
print(f"Length     = {length}")
print(f"Similarity = {similarity:.0f} %")
print(f"Traceback Matrix Size: M[{len(S1)},{len(S2)}]")