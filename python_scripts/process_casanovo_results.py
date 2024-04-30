"""
Script: process_casanovo_results.py
Author: Michael Riffle <mriffle@uw.edu>
Date: April 23, 2024
Description: This script reads one or more tab-delimited text files specified on the command line,
             finds the best PSM line for each distinct peptide based on the highest "search_engine_score[1]",
             and outputs the results to standard out
Usage: python3 process_casanovo_results.py file1.txt file2.txt file3.txt
"""

import sys
import csv
import re

def calculate_error_ppm(expected_mass, observed_mass):
    error_ppm = (observed_mass - expected_mass) / expected_mass * 1000000
    return error_ppm

def process_files(file_paths):
    peptide_data = {}
    peptide_counts = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            # Skip lines until the header line starting with 'PSH' is found
            for row in reader:
                if row[0].startswith('PSH'):
                    headers = row
                    break
            else:
                print(f"Error: No header line starting with 'PSH' found in file: {file_path}")
                return

            # Cache the column indices
            try:
                sequence_index = headers.index('sequence')
                charge_index = headers.index('charge')
                score_index = headers.index('search_engine_score[1]')
                expected_mass_index = headers.index('calc_mass_to_charge')
                observed_mass_index = headers.index('exp_mass_to_charge')
            except ValueError as e:
                print(f"Error: Missing expected column in file: {file_path}")
                print(f"Column not found: {str(e)}")
                return

            # Process each row of data
            for row in reader:
                if not row[0].startswith('PSM'):
                    continue

                sequence = re.sub(r'[^A-Z]', '', row[sequence_index])
                charge = row[charge_index]
                score = float(row[score_index])
                expected_mass = float(row[expected_mass_index])
                observed_mass = float(row[observed_mass_index])
                mz_ppm_error = calculate_error_ppm(expected_mass, observed_mass)

                # Count the number of occurrences of each peptide sequence
                peptide_counts[sequence] = peptide_counts.get(sequence, 0) + 1

                if sequence not in peptide_data or score > peptide_data[sequence]['score']:
                    peptide_data[sequence] = {
                        'charge': charge,
                        'score': score,
                        'file': file_path,
                        'mz_ppm_error': mz_ppm_error,
                        'num_spectra': peptide_counts[sequence]
                    }
                else:
                    peptide_data[sequence]['num_spectra'] = peptide_counts[sequence]

    # Output the results
    print("peptide_sequence\tcharge\tsearch_engine_score[1]\tfile\tmz_ppm_error\tnum_spectra")
    for peptide, data in peptide_data.items():
        print(f"{peptide}\t{data['charge']}\t{data['score']}\t{data['file']}\t{data['mz_ppm_error']:.2f}\t{data['num_spectra']}")

def main():
    """
    Main function to process command-line arguments and call the process_files function.
    """
    # Check if file paths are provided as command-line arguments
    if len(sys.argv) < 2:
        print("Please provide one or more file paths as command-line arguments.")
        return

    # Process the files
    process_files(sys.argv[1:])

if __name__ == '__main__':
    main()