# We are inside a courtroom, Mihai is accused of plagiarism.
# The lawyer must proove his innocence, and asks you for your help
# as a specialist in infractions. You take 2 poetries: one of M. Eminescu
# and the other fo N. Stanescu. You use these 2 texts to create transition
# matrices, which is able to capture  the probability of transition
# between word. You combine the matrices into a Log LikelyHood Matrix.
# We use the LLM to scan, by using a sliding window, the text of Mihai, the 
# guy accused of plagiarism.
# NEgative scores -> one model
# POsitive scores -> other model
# Zero values -> Neither
# Simulate the text of Mihai using ChatGPT using AI and ask it to make it a
# combination between 2 poetries of Eminescu & Stanescu
#At the end of this process you should be able to tell which parts of the text
# in written by Eminescu or by Stanescu or by neither

import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import re

def tokenize(text):
    words = re.findall(r'\b\w+\b', text.lower())
    return words

def build_transition_matrix(words):
    transitions = defaultdict(int)
    totals = defaultdict(int)
    for i in range(len(words)-1):
        a, b = words[i], words[i+1]
        transitions[(a,b)] += 1
        totals[a] += 1
    vocab = sorted(set(words))
    matrix = pd.DataFrame(index=vocab, columns=vocab, data=0.0)
    for (a,b), count in transitions.items():
        matrix.loc[a,b] = count / totals[a]
    return matrix.fillna(0)

poem_eminescu = """
Mai am un singur dor să mă lăsați să mor
La marginea mării-ntr-un tainic loc
Sub un salcâm cu flori
Să mi se-nchidă ochii cu-al cerului noroc
"""

poem_stanescu = """
Se apropie seara cu pași de lumină
Cuvintele tac și gândul alunecă lin
Între real și vis o frântură de vină
Se naște poezia din abisul divin
"""

mihai_text = """
Să mă lăsați să mor cu pași de lumină
Sub un salcâm cu flori, gândul alunecă lin
Cuvintele tac, și mi se-nchid ochii
Se naște dorul din abisul divin
"""

words_eminescu = tokenize(poem_eminescu)
words_stanescu = tokenize(poem_stanescu)
words_mihai = tokenize(mihai_text)

matrix_eminescu = build_transition_matrix(words_eminescu)
matrix_stanescu = build_transition_matrix(words_stanescu)

vocab = sorted(set(matrix_eminescu.index).union(set(matrix_stanescu.index)))
matrix_eminescu = matrix_eminescu.reindex(index=vocab, columns=vocab, fill_value=0)
matrix_stanescu = matrix_stanescu.reindex(index=vocab, columns=vocab, fill_value=0)

epsilon = 1e-10
LLM = np.log2((matrix_stanescu + epsilon) / (matrix_eminescu + epsilon))

def score_window(words, llm, window_size=3):
    results = []
    for i in range(len(words)-window_size):
        score = 0
        for j in range(window_size - 1):
            a, b = words[i + j], words[i + j + 1]
            if a in llm.index and b in llm.columns:
                score += llm.loc[a, b]
        results.append((i, " ".join(words[i:i+window_size]), score))
    return results

results = score_window(words_mihai, LLM)

for idx, segment, score in results:
    label = "Stanescu" if score > 0 else "Eminescu" if score < 0 else "Original Mihai"
    print(f"{segment:40} | Score: {score:6.2f} → {label}")

plt.figure(figsize=(10,6))
sns.heatmap(LLM, cmap="coolwarm", center=0, xticklabels=True, yticklabels=True)
plt.title("Log-Likelihood Matrix between Stanescu / Eminescu")
plt.show()


positions = [i for i, _, _ in results]
scores = [s for _, _, s in results]

plt.figure(figsize=(12, 4))
plt.plot(positions, scores, marker='o')
plt.axhline(0, color='gray', linestyle='--')
plt.title("Sliding Window Log-Likelihood Score on Mihai's Text")
plt.xlabel("Window Position")
plt.ylabel("LLM Score")
plt.grid(True)
plt.show()
