import pandas as pd

def load_matrix(file_path):
    return pd.read_csv(file_path, header=None).values

def load_vector(file_path):
    return pd.read_csv(file_path, header=None).values.flatten()

def predict_steps(matrix, vector, steps=5):
    results = [vector.copy()]
    for _ in range(steps):
        vector = matrix @ vector
        results.append(vector.copy())
    return results

def main():
    matrix = load_matrix('matrix.csv')
    vector = load_vector('vector.csv')

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square.")
    if matrix.shape[1] != vector.shape[0]:
        raise ValueError("Matrix and vector dimensions do not align.")

    results = predict_steps(matrix, vector, steps=5)

    for i, vec in enumerate(results):
        print(f"Step {i}: {vec}")

    pd.DataFrame(results).to_csv("predictions.csv", index=False, header=False)

if __name__ == "__main__":
    main()
