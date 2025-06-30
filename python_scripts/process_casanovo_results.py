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

DELTA_MASS_13C = 1.003355

def calculate_error_ppm(expected_mz, observed_mz, charge):
    min_error_ppm = float('inf')
    for num_13c in range(4):
        adjusted_expected_mz = expected_mz + (num_13c * DELTA_MASS_13C) / charge
        error_ppm = (observed_mz - adjusted_expected_mz) / adjusted_expected_mz * 1000000
        if abs(error_ppm) < abs(min_error_ppm):
            min_error_ppm = error_ppm
    return min_error_ppm


def add_rank_score_to_peptide_data(peptide_data):
    total_peptides = len(peptide_data)
    
    # Extract scores and sort them in descending order (highest score = best rank)
    scores = [data['score'] for data in peptide_data.values()]
    unique_scores = sorted(set(scores), reverse=True)
    
    # Create a mapping from score to rank
    score_to_rank = {}
    current_rank = 1
    
    for score in unique_scores:
        score_to_rank[score] = current_rank

        # Count how many peptides have this score to determine next rank
        count_with_score = scores.count(score)
        current_rank += count_with_score
    
    # Add rank_score to each peptide
    for sequence, data in peptide_data.items():
        rank = score_to_rank[data['score']]
        data['rank_score'] = rank / total_peptides


def process_files(file_paths):
    peptide_data = {}
    peptide_counts = {}
    peptide_peptidoforms = {}

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
                expected_mz_index = headers.index('calc_mass_to_charge')
                observed_mz_index = headers.index('exp_mass_to_charge')
            except ValueError as e:
                print(f"Error: Missing expected column in file: {file_path}")
                print(f"Column not found: {str(e)}")
                return

            # Process each row of data
            for row in reader:
                if not row[0].startswith('PSM'):
                    continue

                peptidoform_sequence = row[sequence_index]
                sequence = re.sub(r'[^A-Z]', '', row[sequence_index])
                charge = int(float(row[charge_index]))
                score = float(row[score_index])
                expected_mz = float(row[expected_mz_index])
                observed_mz = float(row[observed_mz_index])
                mz_ppm_error = calculate_error_ppm(expected_mz, observed_mz, charge)

                # We effectively disable Casanovo’s precursor m/z filtering by adding 1
                # to any negative Casanovo score, thereby ensuring that all scores are
                # in the range [0,1].
                if score < 0:
                    score += 1

                # Count the number of occurrences of each peptide sequence
                peptide_counts[sequence] = peptide_counts.get(sequence, 0) + 1

                # add this peptidoform to the set of peptidoforms for this peptide
                if sequence not in peptide_peptidoforms:
                    peptide_peptidoforms[sequence] = set()
                
                peptide_peptidoforms[sequence].add(peptidoform_sequence)

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

    # add a rank score to each peptide
    add_rank_score_to_peptide_data(peptide_data)

    # Output the results
    print("peptide_sequence\tcharge\tsearch_engine_score[1]\tfile\tmz_ppm_error\tnum_spectra\trank_score\tnum_peptidoforms")
    for peptide, data in peptide_data.items():
        print(f"{peptide}\t{data['charge']}\t{data['score']}\t{data['file']}\t{data['mz_ppm_error']:.2f}\t{data['num_spectra']}\t{data['rank_score']}\t{len(peptide_peptidoforms[peptide])}")

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