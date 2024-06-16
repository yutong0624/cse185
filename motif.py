import argparse
from scipy.stats import binomtest
import numpy as np
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import re
import matplotlib.pyplot as plt
def parse_fasta(file_path):
    """Reads a FASTA file and returns a list of sequences."""
    sequences = []
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence.replace('\n', '').replace(' ', ''))
                sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence.replace('\n', '').replace(' ', ''))
    return sequences



def find_motifs(sequence, motif_length):
    """Finds all motifs of a given length in a sequence."""
    motifs = Counter()

    for i in range(len(sequence) - motif_length + 1):
        motif = sequence[i:i + motif_length]
        motifs[motif] += 1
    return motifs

def count_motifs(sequences, motif_length):
    """Counts motifs across multiple sequences using parallel processing."""
    motif_counts = Counter()
    with ProcessPoolExecutor() as executor:
        results = executor.map(find_motifs, sequences, [motif_length] * len(sequences))
        for result in results:
            motif_counts.update(result)
    return motif_counts

def filter_significant_motifs(motif_counts, background_counts, total_sequences, motif_length, sequences, background_sequences):
    """Filters motifs based on significance using binomial test and background frequencies."""
    significant_motifs = []
    sequence_length = len(sequences[0])
    total_positions = total_sequences * (sequence_length - motif_length + 1)
    background_sequence_length = len(background_sequences[0])
    
    # Determine the cutoff for the top 20% motifs by count
    counts = np.array(list(motif_counts.values()))
    
    for motif, count in motif_counts.items():
        if count < 10:
            continue
        background_motif_count = background_counts.get(motif, 0)
        expected_count = background_motif_count * total_positions / (len(background_sequences) * (background_sequence_length - motif_length + 1))
        # Calculate the expected probability
        expected_prob = expected_count / total_positions
        if expected_prob > 1:
            expected_prob = 1
        p_value = binomtest(count, total_positions, expected_prob, alternative='greater').pvalue
        if p_value < 0.05:
            significant_motifs.append((motif, count, p_value))
    # Sort motifs by p-value in ascending order
    significant_motifs.sort(key=lambda x: x[2])
    return significant_motifs[:5]  # Return top 5 motifs with the lowest p-values

def plot_motif_frequencies(significant_motifs):
    """Plots the frequencies of the top significant motifs."""
    motifs, counts, p_values = zip(*significant_motifs)
    
    plt.figure(figsize=(10, 5))
    plt.bar(motifs, counts, color='blue', edgecolor='black')
    plt.xlabel('Motifs')
    plt.ylabel('Counts')
    plt.title('Top Significant Motifs')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Find and count significant DNA motifs in sequences from a text file compared to a background model.")
    parser.add_argument('target_file', type=str, help='Input text file with target sequences.')
    parser.add_argument('background_file', type=str, help='Input FASTA file with background sequences.')
    parser.add_argument('-l', '--length', type=int, required=True, help='Length of the motifs to find')
    args = parser.parse_args()

    target_sequences = parse_fasta(args.target_file)
   
    background_sequences = parse_fasta(args.background_file)



    target_counts = count_motifs(target_sequences, args.length)
    background_counts = count_motifs(background_sequences, args.length)

    print(f"Target counts: {len(target_counts)} motifs")
    print(f"Background counts: {len(background_counts)} motifs")

    significant_motifs = filter_significant_motifs(target_counts, background_counts, len(target_sequences), args.length, target_sequences, background_sequences)

    if not significant_motifs:
        print("No significant motifs found.")
    else:
        for motif, count, p_value in significant_motifs:
            print(f"Motif: {motif}, Count: {count}, p-value: {p_value}")
        plot_motif_frequencies(significant_motifs)

if __name__ == "__main__":
    main()
