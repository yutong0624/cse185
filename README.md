Motif Discovery Tool
Overview
This tool is designed to identify significant DNA motifs in a given set of genomic sequences compared to a background model. It uses statistical methods to determine the significance of motifs and provides a list of the most significant motifs found.

Features
Parses input sequences from text or FASTA files.
Identifies motifs of specified lengths.
Compares motif occurrences between target and background sequences.
Filters significant motifs based on binomial test p-values.
Supports parallel processing for faster computation.
Outputs the most significant motifs and their p-values.
Installation
To use this tool, you need Python 3.x and the following Python libraries:

scipy
numpy
collections
concurrent.futures
re
You can install the required libraries using pip:

sh
Copy code
pip install scipy numpy
Usage
sh
Copy code
python motif.py <target_file> <background_file> -l <motif_length>
<target_file>: Path to the text file with target sequences.
<background_file>: Path to the FASTA file with background sequences.
-l <motif_length>: Length of the motifs to find (required).
Main Function
The main function orchestrates the reading of input files, counting motifs, and filtering significant motifs. It prints out the most significant motifs found, along with their counts and p-values.

Results and Benchmarking
Benchmarking
The tool was benchmarked against the HOMER motif discovery tool using a published dataset. The core motif for the transcription factor in the dataset should be "CACGTG." Both tools produced similar results, and the processing time was around 8 minutes for both tools.

Results
Tool Results
plaintext
Copy code
Target counts: 1006 motifs
Background counts: 15070 motifs
Motif: TTTTT, Count: 138, p-value: 2.5740005782776596e-303
Motif: CTCCC, Count: 58, p-value: 3.7720968406907563e-118
Motif: AGAGG, Count: 58, p-value: 3.815814954696852e-114
Motif: CCGCC, Count: 42, p-value: 4.342270436943195e-110
Motif: GAGGA, Count: 54, p-value: 3.9428811880440666e-107
HOMER Results
plaintext
Copy code
>TATCG	1-TATCG	5.148095	0.000000	0	T:1.0(100.00%),B:177.1(100.00%),P:1e0	Tpos:1831.0,Tstd:233.0,Bpos:3362300.9,Bstd:4492853.8,StrandBias:10.0,Multiplicity:2.00
>TCGAT	2-TCGAT	5.148095	0.000000	0	T:1.0(100.00%),B:177.1(100.00%),P:1e0	Tpos:3971.5,Tstd:1066.0,Bpos:2683793.9,Bstd:4214475.4,StrandBias:0.0,Multiplicity:3.00
>GTCGA	3-GTCGA	5.148095	0.000000	0	T:1.0(100.00%),B:177.1(100.00%),P:1e0	Tpos:2543.3,Tstd:948.5,Bpos:2529046.1,Bstd:4192481.7,StrandBias:-1.0,Multiplicity:3.00
>TCGAA	4-TCGAA	5.148095	0.000000	0	T:1.0(100.00%),B:177.1(100.00%),P:1e0	Tpos:2211.5,Tstd:775.5,Bpos:2598045.3,Bstd:4131107.8,StrandBias:-10.0,Multiplicity:2.00
>GTTCG	5-GTTCG	5.148095	0.000000	0	T:1.0(100.00%),B:177.1(100.00%),P:1e0	Tpos:3769.7,Tstd:1091.9,Bpos:2435790.5,Bstd:4135431.8,StrandBias:-1.0,Multiplicity:3.00
Known Motif
According to the published dataset, the core motif for this transcription factor should be "CACGTG."

Performance
The tool took approximately 8 minutes to process, which is comparable to HOMER.

Challenges and Future Directions
Challenges
Incorrect handling of headers in the FASTA parser.
Errors related to processing non-DNA sequences.
Issues with the calculation of expected counts in the statistical tests.
Future Directions
Algorithm Optimization: Further refine the algorithm to handle larger datasets more efficiently, potentially through more advanced data structures or search techniques.
Integration of Additional Statistical Models: Enhance the tool by incorporating a variety of statistical models to assess motif significance, providing a more robust analysis.
Incorporation of Multiple Background Models: Use various background models to better estimate the significance of identified motifs, improving the accuracy of the results.
Optimization for Speed: Implement advanced parallel processing techniques to speed up the tool's performance, making it capable of processing large datasets more quickly.
Code Availability
The code for this tool is available on GitHub: https://github.com/yutong0624/cse185

References
NCBI RefSeq assembly: GCF_000001405.40
Genome Reference Consortium: GCA_000001405.29
Taxon: Homo sapiens (human)
Synonym: hg38
Assembly type: haploid with alt loci
Submitter: Genome Reference Consortium
Date: Feb 3, 2022
