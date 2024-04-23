"""
Script: process_comet_results.py
Author: Michael Riffle <mriffle@uw.edu>
Date: April 23, 2024
Description: This script reads one or more tab-delimited text files specified on the command line,
             finds the best line for each distinct peptide based on the lowest "e-value", and outputs
             the results to standard out.
Usage: python3 process_comet_results.py file1.txt file2.txt file3.txt
"""

import sys
import csv

def process_files(file_paths):

    peptide_data = {}

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            
            # Skip the first line
            next(file)
            
            reader = csv.reader(file, delimiter='\t')
            
            # Get the headers and cache the column indices
            headers = next(reader)

            try:
                plain_peptide_index = headers.index('plain_peptide')
                charge_index = headers.index('charge')
                e_value_index = headers.index('e-value')
                protein_index = headers.index('protein')
            except ValueError as e:
                print(f"Error: Missing expected column in file: {file_path}")
                print(f"Column not found: {str(e)}")
                return
            
            # Process each row of data
            for row in reader:
                plain_peptide = row[plain_peptide_index]
                charge = row[charge_index]
                e_value = float(row[e_value_index])
                protein = row[protein_index]
                
                if plain_peptide not in peptide_data or e_value < peptide_data[plain_peptide]['e_value']:
                    peptide_data[plain_peptide] = {
                        'charge': charge,
                        'e_value': e_value,
                        'protein': protein,
                        'file': file_path
                    }

    # Output the results
    print("plain_peptide\tcharge\te-value\tprotein\tfile")
    for peptide, data in peptide_data.items():
        print(f"{peptide}\t{data['charge']}\t{data['e_value']}\t{data['protein']}\t{data['file']}")

def main():
    # Check if file paths are provided as command-line arguments
    if len(sys.argv) < 2:
        print("Please provide one or more file paths as command-line arguments.")
        return
    
    # Process the files
    process_files(sys.argv[1:])

if __name__ == '__main__':
    main()
