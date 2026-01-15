# Use the trasnsition mstrix from the JSON file in order to
# synthesise new sequences of text based on the transition matrix

import json
import random

def load_transition_matrix(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def generate_sequence_from_matrix(matrix, start_token, length=30):
    sequence = [start_token]
    current = start_token

    for _ in range(length - 1):
        if current not in matrix:
            break
        next_tokens = list(matrix[current].keys())
        probabilities = list(matrix[current].values())
        current = random.choices(next_tokens, weights=probabilities)[0]
        sequence.append(current)

    return sequence

def main():
    # Text
    matrix_file_text = 'word_transition_matrix.json'
    start_text = 'the'
    matrix_text = load_transition_matrix(matrix_file_text)
    sequence_text = generate_sequence_from_matrix(matrix_text, start_token=start_text, length=50)

    # DNA
    matrix_file_dna = 'dna_transition_matrix.json'
    start_dna = 'A'
    matrix_dna = load_transition_matrix(matrix_file_dna)
    sequence_dna = generate_sequence_from_matrix(matrix_dna, start_token=start_dna, length=50)

    # Print
    print("\nGenerated sequence (TEXT):")
    print(' '.join(sequence_text))

    print("\nGenerated sequence (DNA):")
    print(''.join(sequence_dna))

if __name__ == '__main__':
    main()



