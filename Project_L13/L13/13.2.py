import random
import json
import re
from collections import defaultdict

def save_matrix_to_json(matrix, filename):
    with open(filename, 'w') as f:
        json.dump({k: dict(v) for k, v in matrix.items()}, f, indent=4)

def load_transition_matrix(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def generate_random_dna(length=50):
    return ''.join(random.choices('ACGT', k=length))

def compute_transition_matrix_dna(sequence):
    matrix = defaultdict(lambda: defaultdict(int))
    total_transitions = defaultdict(int)

    for i in range(len(sequence) - 1):
        current = sequence[i]
        next_ = sequence[i + 1]
        matrix[current][next_] += 1
        total_transitions[current] += 1

    for curr in matrix:
        for next_ in matrix[curr]:
            matrix[curr][next_] /= total_transitions[curr]

    return matrix

def generate_sequence_dna(matrix, start='A', length=50):
    sequence = [start]
    current = start

    for _ in range(length - 1):
        if current not in matrix:
            break
        next_chars = list(matrix[current].keys())
        probabilities = list(matrix[current].values())
        current = random.choices(next_chars, weights=probabilities)[0]
        sequence.append(current)

    return ''.join(sequence)

def generate_random_english_text():
    sample_text = """
    The quick brown fox jumps over the lazy dog. A small cat watches silently. 
    Then suddenly, the wind blows and everything changes. Nature whispers, but 
    we rarely listen. Curiosity sparks discovery and discovery leads to progress.
    """
    return sample_text[:300]

def tokenize_words(text):
    cleaned = re.sub(r'[^\w\s\.\!\?]', '', text.lower())
    words = cleaned.split()
    return words

def compute_transition_matrix_words(words):
    matrix = defaultdict(lambda: defaultdict(int))
    total_transitions = defaultdict(int)

    for i in range(len(words) - 1):
        current = words[i]
        next_ = words[i + 1]
        matrix[current][next_] += 1
        total_transitions[current] += 1

    for curr in matrix:
        for next_ in matrix[curr]:
            matrix[curr][next_] /= total_transitions[curr]

    return matrix

def generate_sequence_text(matrix, start='the', length=20):
    sequence = [start]
    current = start

    for _ in range(length - 1):
        if current not in matrix:
            break
        next_words = list(matrix[current].keys())
        probabilities = list(matrix[current].values())
        current = random.choices(next_words, weights=probabilities)[0]
        sequence.append(current)

    return ' '.join(sequence)

def main():
    dna_seq = generate_random_dna()
    dna_matrix = compute_transition_matrix_dna(dna_seq)
    save_matrix_to_json(dna_matrix, 'dna_transition_matrix.json')
    print("DNA Transition Matrix saved to dna_transition_matrix.json")
    print("Original DNA Sequence:", dna_seq)

    # Generate new DNA sequence
    generated_dna = generate_sequence_dna(dna_matrix, start=dna_seq[0], length=50)
    print("Generated DNA Sequence:", generated_dna)

    # --- TEXT ---
    text = generate_random_english_text()
    words = tokenize_words(text)
    word_matrix = compute_transition_matrix_words(words)
    save_matrix_to_json(word_matrix, 'word_transition_matrix.json')
    print("Word Transition Matrix saved to word_transition_matrix.json")

    # Generate new text sequence
    generated_text = generate_sequence_text(word_matrix, start='the', length=30)
    print("Generated Text:", generated_text)

if __name__ == '__main__':
    main()
