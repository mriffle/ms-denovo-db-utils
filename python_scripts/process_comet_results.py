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

def is_n_tryptic(modified_peptide):
    return 1 if modified_peptide.startswith(('R', 'K')) else 0

def is_c_tryptic(plain_peptide):
    return 1 if plain_peptide.endswith(('R', 'K')) else 0

def calculate_error_ppm(expected_mass, observed_mass):
    error_ppm = (observed_mass - expected_mass) / expected_mass * 1000000
    return error_ppm

def process_files(file_paths):
    peptide_data = {}
    peptide_counts = {}
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
                modified_peptide_index = headers.index('modified_peptide')
                calc_neutral_mass_index = headers.index('calc_neutral_mass')
                exp_neutral_mass_index = headers.index('exp_neutral_mass')
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
                modified_peptide = row[modified_peptide_index]
                calc_neutral_mass = float(row[calc_neutral_mass_index])
                exp_neutral_mass = float(row[exp_neutral_mass_index])
                
                # Increment the count for the plain_peptide
                peptide_counts[plain_peptide] = peptide_counts.get(plain_peptide, 0) + 1
                
                if plain_peptide not in peptide_data or e_value < peptide_data[plain_peptide]['e_value']:
                    peptide_data[plain_peptide] = {
                        'charge': charge,
                        'e_value': e_value,
                        'protein': protein,
                        'file': file_path,
                        'tryptic_n': is_n_tryptic(modified_peptide),
                        'tryptic_c': is_c_tryptic(plain_peptide),
                        'ppm_error': calculate_error_ppm(calc_neutral_mass, exp_neutral_mass)
                    }
    
    # Output the results
    print("plain_peptide\tcharge\te-value\tprotein\tfile\ttryptic_n\ttryptic_c\tnum_spectra\tppm_error")
    for peptide, data in peptide_data.items():
        num_spectra = peptide_counts[peptide]
        print(f"{peptide}\t{data['charge']}\t{data['e_value']}\t{data['protein']}\t{data['file']}\t{data['tryptic_n']}\t{data['tryptic_c']}\t{num_spectra}\t{data['ppm_error']:.2f}")

def main():
    # Check if file paths are provided as command-line arguments
    if len(sys.argv) < 2:
        print("Please provide one or more file paths as command-line arguments.")
        return
    
    # Process the files
    process_files(sys.argv[1:])

if __name__ == '__main__':
    main()