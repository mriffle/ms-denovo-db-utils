"""
Script: collate_into_fasta.py
Author: Michael Riffle <mriffle@uw.edu>
Date: April 23, 2024
Description: This script reads one or more tab-delimited text files specified on the command line,
             takes the peptide sequence in the first column each file, collects all distinct
             peptides and oupts the final list in FASTA format using:
               >peptide_sequence
               peptide_sequence
Usage: python3 collate_into_fasta.py file1.txt file2.txt file3.txt
"""

import sys

def process_files(file_paths):
    peptide_sequences = set()

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            # Skip the header line
            next(file)

            for line in file:
                columns = line.strip().split('\t')
                if columns:
                    peptide_sequence = columns[0]
                    peptide_sequences.add(peptide_sequence)

    return peptide_sequences

def write_fasta_output(peptide_sequences):
    for peptide_sequence in peptide_sequences:
        print(f">{peptide_sequence}")
        print(peptide_sequence)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide at least one file path as a command line argument.")
        sys.exit(1)

    file_paths = sys.argv[1:]
    peptide_sequences = process_files(file_paths)
    write_fasta_output(peptide_sequences)