#Design an application using the AI which contains GUI which allows the user 
#select a FASTA file. The content of the file should be analyzed by using a 
#sliding window of 30 positions. The content of each sliding window should be 
#used in order to extract the relative frequencies of the symbols found in the 
#alphabet of the sequence. Thus your input should be the DNA seq from the 
#FASTA file, and the output should be the values of relative frequencies of 
#each symbol in the alphabet of the sequence. Translate in lines on a chart, 
#thus your chart in the case of DNA should have 4 lines which reflect the 
#values found over the sequence

from tkinter import filedialog
from tkinter import Tk, Label, Button, messagebox
import matplotlib.pyplot as plt

def relative_frequencies(sequence, window_size):
    bases = ['A', 'T', 'C', 'G']
    freq = {base: [] for base in bases}

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        for base in bases:
            freq[base].append(round(window.count(base) / window_size, 4))

    return freq

def show_chart(freq):
    for base in freq:
        plt.plot(freq[base], label=base)
    plt.title("Relative Frequencies of A, T, C, G (Window = 30)")
    plt.xlabel("Window Start Position")
    plt.ylabel("Relative Frequency")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def open_file():
    filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.txt")])
    if not filepath:
        return
    try:
        with open(filepath, "r") as file:
            lines = file.readlines()
            sequence = "".join([line.strip().upper() for line in lines if not line.startswith(">")])
            sequence = "".join([base for base in sequence if base in ['A', 'T', 'C', 'G']])

            if len(sequence) < 30:
                messagebox.showerror("Error", "Sequence too short.")
                return

            freq = relative_frequencies(sequence, 30)
            show_chart(freq)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file:\n{e}")

root = Tk()
root.title("Sliding Window A/T/C/G")
root.geometry("400x200")

label = Label(root, text="Select a FASTA file to analyze", font=("Arial", 12))
label.pack(pady=20)

button = Button(root, text="Open FASTA File", font=("Arial", 12), command=open_file)
button.pack()

root.mainloop()
