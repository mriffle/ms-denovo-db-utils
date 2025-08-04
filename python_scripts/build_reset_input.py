import os
import sys
import math

def read_peptide_file(file_path):
    peptide_map = {}

    with open(file_path, 'r') as file:
        header = file.readline().strip().split('\t')

        for line in file:
            columns = line.strip().split('\t')
            plain_peptide = columns[0]
            peptide_data = {header[i]: columns[i] for i in range(1, len(header))}
            peptide_map[plain_peptide] = peptide_data

    return peptide_map

def read_diamond_file(file_path, library_decoy_prefix):
    """Read Diamond results in outfmt 6 format"""
    column_headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    peptide_map = {}
    peptide_protein_map = {}

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            # Skip empty lines and comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            columns = line.split('\t')

            if len(columns) != len(column_headers):
                raise ValueError(f"Line {line_num}: Expected {len(column_headers)} columns, but found {len(columns)}")

            # In Diamond format, qseqid is the query (peptide sequence)
            peptide_sequence = columns[0]
            protein_name = columns[1]

            # Create data dictionary with all the alignment information
            peptide_data = {column_headers[i]: columns[i] for i in range(len(column_headers))}
            
            # If we have multiple hits for the same peptide, keep the best one (lowest e-value)
            if peptide_sequence in peptide_map:
                existing_evalue = float(peptide_map[peptide_sequence]['evalue'])
                new_evalue = float(peptide_data['evalue'])
                if new_evalue < existing_evalue:
                    peptide_map[peptide_sequence] = peptide_data
                    peptide_protein_map[peptide_sequence] = set()
            else:
                peptide_map[peptide_sequence] = peptide_data
                peptide_protein_map[peptide_sequence] = set()

            peptide_protein_map[peptide_sequence].add(protein_name)
        
    # remove peptides that map to both decoys and non decoys in the annotated db search
    peptides_to_remove = set()
    for peptide in peptide_protein_map:
        any_starts_with_decoy = any(s.startswith(library_decoy_prefix) for s in peptide_protein_map[peptide])
        any_not_starts_with_decoy = any(not s.startswith(library_decoy_prefix) for s in peptide_protein_map[peptide])

        if any_starts_with_decoy and any_not_starts_with_decoy:
            peptides_to_remove.add(peptide)

    filtered_peptide_map = {k: v for k, v in peptide_map.items() if k not in peptides_to_remove}

    return filtered_peptide_map

import gzip

def augment_peptide_map(peptide_map, fasta_file_path):
    """Extract subject sequences from FASTA file based on Diamond alignment coordinates"""
    protein_peptide_map = {}
    for peptide, data in peptide_map.items():
        protein_name = data['sseqid']
        if protein_name not in protein_peptide_map:
            protein_peptide_map[protein_name] = []
        protein_peptide_map[protein_name].append(peptide)

    # Check if file is gzipped and open accordingly
    if fasta_file_path.endswith('.gz'):
        fasta_file = gzip.open(fasta_file_path, 'rt')
    else:
        fasta_file = open(fasta_file_path, 'r')
    
    try:
        protein_name = None
        protein_sequence = ''

        for line in fasta_file:
            line = line.strip()

            if line.startswith('>'):
                if protein_name and protein_name in protein_peptide_map:
                    for peptide in protein_peptide_map[protein_name]:
                        start = int(peptide_map[peptide]['sstart']) - 1
                        end = int(peptide_map[peptide]['send'])
                        subject_sequence = protein_sequence[start:end]
                        peptide_map[peptide]['ssequence'] = subject_sequence

                protein_name = line[1:].split()[0]
                protein_sequence = ''
            else:
                protein_sequence += line

        if protein_name and protein_name in protein_peptide_map:
            for peptide in protein_peptide_map[protein_name]:
                start = int(peptide_map[peptide]['sstart']) - 1
                end = int(peptide_map[peptide]['send'])
                subject_sequence = protein_sequence[start:end]
                peptide_map[peptide]['ssequence'] = subject_sequence
    
    finally:
        fasta_file.close()

    for peptide, data in peptide_map.items():
        if 'ssequence' not in data:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")

    return peptide_map

def build_diamond_hit_query_sequence_map(peptides, diamond_map):
    """
    Build a dict mapping library_peptide_sequence to the set of
    casanovo/comet peptides that mapped to that library hit.

    Parameters:
        peptides (iterable): List or iterable of peptide strings to process.
        diamond_map (dict): Mapping from peptide => diamond_data dict.

    Returns:
        dict: {ssequence => set of peptides}
    """
    diamond_hit_query_sequence_map = {}

    for peptide in peptides:
        if peptide not in diamond_map:
            continue

        diamond_data = diamond_map[peptide]

        ssequence = diamond_data.get('ssequence')
        if ssequence is None:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")
        
        if ssequence not in diamond_hit_query_sequence_map:
            diamond_hit_query_sequence_map[ssequence] = set()
        
        diamond_hit_query_sequence_map[ssequence].add(peptide)

    return diamond_hit_query_sequence_map

def output_peptide_data_for_reset(comet_map, casanovo_map, diamond_map, library_decoy_prefix, comet_decoy_prefix):
    peptides = set(comet_map.keys()) | set(casanovo_map.keys())
    missing_peptides = peptides - set(diamond_map.keys())

    if missing_peptides:
        warning_message = f"Warning: The following {len(missing_peptides)} peptides are missing from diamond_map and will be skipped: " + ", ".join(missing_peptides)
        print(warning_message, file=sys.stderr)

    # create a dict of { library_peptide_sequence => set(casanovo/comet
    # peptides that mapped to that library hit) }
    diamond_hit_query_sequence_map = build_diamond_hit_query_sequence_map(peptides, diamond_map)

    header = [
        "SpecId", "Label", "ScanNr",
        "database_peptide_length", "max_diamond_bitscore", "max_diamond_perc_identity", "num_casanovo_peptides", "num_comet_peptides",
        "casanovo_num_spectra", "casanovo_best_score", "casanovo_ppm_error", "casanovo_num_peptidoforms",
        "comet_num_spectra", "comet_n_tryptic", "comet_c_tryptic", "comet_best_score", "comet_ppm_error", "comet_num_peptidoforms",
        "combined_rank_score"
    ]

    header.extend(["Peptide", "Proteins"])
    print("\t".join(header))

    scan_nr = 1

    # iterate over annotated diamond hit peptides, each is a row in RESET input
    for library_hit_peptide, casanovo_comet_search_peptides in diamond_hit_query_sequence_map.items():

        # best diamond hit data
        best_diamond_peptide_length = 0
        best_diamond_bit_score = 0
        best_diamond_perc_identity = 0
        best_diamond_label = 1
        best_diamond_ssequence = ''
        best_diamond_protein = ''

        # best casanovo hit data
        casanovo_num_spectra = 0
        casanovo_num_peptides = 0
        casanovo_best_score = 0
        casanovo_ppm_error = 0
        casanovo_num_peptidoforms = 0
        casanovo_best_rank_score = 2

        # best comet hit data
        comet_num_spectra = 0
        comet_num_peptides = 0
        comet_n_tryptic = 0
        comet_c_tryptic = 0
        comet_best_score = 0
        comet_ppm_error = 0
        comet_num_peptidoforms = 0
        comet_best_rank_score = 2
        comet_best_is_decoy = False

        for peptide in casanovo_comet_search_peptides:
            diamond_data = diamond_map[peptide]
            casanovo_data = casanovo_map.get(peptide, {})
            comet_data = comet_map.get(peptide, {})
        
            # check for best diamond hit
            diamond_bitscore = float(diamond_data.get('bitscore', 0))
            if diamond_bitscore > best_diamond_bit_score:
                database_peptide_length = int(diamond_data.get('send', 0)) - int(diamond_data.get('sstart', 0)) + 1

                best_diamond_bit_score = diamond_bitscore
                best_diamond_peptide_length = database_peptide_length
                best_diamond_perc_identity = diamond_data.get('pident', 0)
                best_diamond_label = -1 if diamond_data.get('sseqid', '').startswith(library_decoy_prefix) else 1
                best_diamond_ssequence = diamond_data.get('ssequence')
                best_diamond_protein = diamond_data.get('sseqid', '')

            # add in casanovo data
            if len(casanovo_data) > 0:

                casanovo_num_spectra += int(casanovo_data.get('num_spectra'))
                casanovo_num_peptidoforms += int(casanovo_data.get('num_peptidoforms'))
                casanovo_num_peptides += 1

                # check for best casanovo hit
                casanovo_score = float(casanovo_data.get('search_engine_score[1]'))
                if casanovo_score > casanovo_best_score:
                    casanovo_best_score = casanovo_score
                    casanovo_ppm_error = casanovo_data.get('mz_ppm_error')
                    casanovo_best_rank_score = casanovo_data.get('rank_score')

            # add in comet data
            if len(comet_data) > 0:

                # skip this comet hit if it's a decoy hit hitting a target library hit
                comet_is_decoy = comet_data.get('is_decoy') == '1'

                if comet_is_decoy and best_diamond_label == 1:
                    warning_message = f"Warning: Ignoring decoy comet hit for {peptide} that matched target in library: {best_diamond_protein}"
                    print(warning_message, file=sys.stderr)
                
                else:
                    comet_num_spectra += int(comet_data.get('num_spectra'))
                    comet_num_peptidoforms += int(comet_data.get('num_peptidoforms'))
                    comet_num_peptides += 1

                    # check for best comet hit
                    comet_score = math.log10(1 + (1 / (float(comet_data.get('e-value')) + 1E-20)))
                    if comet_score > comet_best_score:
                        comet_best_score = comet_score
                        comet_ppm_error = comet_data.get('mz_ppm_error')
                        comet_best_rank_score = comet_data.get('rank_score')
                        comet_n_tryptic = comet_data.get('tryptic_n')
                        comet_c_tryptic = comet_data.get('tryptic_c')
                        comet_best_is_decoy = comet_is_decoy
                

        # skip this row if there were no casanovo and no comet hits (after filtering out decoys that matched library targets)
        if casanovo_num_spectra == 0 and comet_num_spectra == 0:
            continue

        # sanity check that the best comet hit is never a decoy matching a annotated db target
        if comet_best_is_decoy and best_diamond_label == 1:
             raise ValueError(f"Found situation where best comet hit is a decoy and best diamond hit is a target, stopping.")

        spec_id = library_hit_peptide + '_1'
        combined_rank_score = str(4 - float(comet_best_rank_score) - float(casanovo_best_rank_score))
        
        row = [
            spec_id, best_diamond_label, scan_nr, best_diamond_peptide_length, best_diamond_bit_score, best_diamond_perc_identity,
            casanovo_num_peptides, comet_num_peptides,
            casanovo_num_spectra, casanovo_best_score, casanovo_ppm_error, casanovo_num_peptidoforms,
            comet_num_spectra, comet_n_tryptic, comet_c_tryptic, comet_best_score, comet_ppm_error, comet_num_peptidoforms,
            combined_rank_score,
            library_hit_peptide, best_diamond_protein
        ]

        print("\t".join(str(value) for value in row))
        scan_nr += 1


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python build_reset_input.py <comet_results_file> <casanovo_results_file> <diamond_results_file> <fasta_file> <annotated db decoy prefix> <comet decoy prefix>")
        sys.exit(1)

    comet_results_file = sys.argv[1]
    casanovo_results_file = sys.argv[2]
    diamond_results_file = sys.argv[3]
    fasta_file = sys.argv[4]
    library_decoy_prefix = sys.argv[5]
    comet_decoy_prefix = sys.argv[6]

    comet_map = read_peptide_file(comet_results_file)
    casanovo_map = read_peptide_file(casanovo_results_file)
    diamond_map = read_diamond_file(diamond_results_file, library_decoy_prefix)
    diamond_map = augment_peptide_map(diamond_map, fasta_file)

    output_peptide_data_for_reset(comet_map, casanovo_map, diamond_map, library_decoy_prefix)