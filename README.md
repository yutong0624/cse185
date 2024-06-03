# cse185
Description: Motif Finder is a Python command-line utility designed to identify and count DNA motifs in genomic sequences. The tool allows users to specify the motif length and provides options for visualizing motif frequencies through histograms.

Motif Identification: Identify all possible motifs of a specified length in given DNA sequences.
Motif Counting: Count the occurrences of each motif across multiple sequences.
Histogram Plotting: Visualize the frequency distribution of identified motifs.

File:
motif.py: The main Python script for motif finding and counting.
queries2.txt: Example input file containing DNA sequences.

commands: python motif.py queries2.txt -l 5 -p
This command will identify and count all motifs of length 5 in the sequences provided in queries2.txt and plot a histogram of their frequencies.
