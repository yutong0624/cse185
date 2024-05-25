import argparse
import numpy as np
import matplotlib.pyplot as plt

def parse_txt(file_path):
    """Reads a text file and returns a list of sequences, one per line."""
    with open(file_path, 'r') as file:
        sequences = [line.strip() for line in file if line.strip()]
    return sequences

def find_motifs(sequence, motif_length):
    """Finds all motifs of a given length in a sequence."""
    motifs = {}
    for i in range(len(sequence) - motif_length + 1):
        motif = sequence[i:i + motif_length]
        if motif in motifs:
            motifs[motif] += 1
        else:
            motifs[motif] = 1
    return motifs

def count_motifs(sequences, motif_length):
    """Counts motifs across multiple sequences."""
    motif_counts = {}
    for seq in sequences:
        motifs = find_motifs(seq, motif_length)
        for motif, count in motifs.items():
            if motif in motif_counts:
                motif_counts[motif] += count
            else:
                motif_counts[motif] = count
    return motif_counts

def plot_histogram(motif_counts):
    """Plots a histogram of motif frequencies."""
    plt.figure(figsize=(10, 5))
    values = list(motif_counts.values())
    plt.hist(values, bins=np.arange(1, max(values) + 1, 1), alpha=0.75, color='blue', edgecolor='black')
    plt.title('Histogram of Motif Frequencies')
    plt.xlabel('Frequency')
    plt.ylabel('Number of Motifs')
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Find and count DNA motifs in sequences from a text file.")
    parser.add_argument('txt_file', type=str, help='Input text file with sequences.')
    parser.add_argument('-l', '--length', type=int, required=True, help='Length of the motifs to find')
    parser.add_argument('-p', '--plot', action='store_true', help='Plot histogram of motif frequencies')
    args = parser.parse_args()

    sequences = parse_txt(args.txt_file)
    motif_counts = count_motifs(sequences, args.length)

    for motif, count in motif_counts.items():
        print(f"Motif: {motif}, Count: {count}")

    if args.plot:
        plot_histogram(motif_counts)

if __name__ == "__main__":
    main()


