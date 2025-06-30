import os
import sys

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

def read_diamond_file(file_path):
    """Read Diamond results in outfmt 6 format"""
    column_headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    peptide_map = {}

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
            
            # Create data dictionary with all the alignment information
            peptide_data = {column_headers[i]: columns[i] for i in range(len(column_headers))}
            
            # If we have multiple hits for the same peptide, keep the best one (lowest e-value)
            if peptide_sequence in peptide_map:
                existing_evalue = float(peptide_map[peptide_sequence]['evalue'])
                new_evalue = float(peptide_data['evalue'])
                if new_evalue < existing_evalue:
                    peptide_map[peptide_sequence] = peptide_data
            else:
                peptide_map[peptide_sequence] = peptide_data

    return peptide_map

def augment_peptide_map(peptide_map, fasta_file_path):
    """Extract subject sequences from FASTA file based on Diamond alignment coordinates"""
    protein_peptide_map = {}
    for peptide, data in peptide_map.items():
        protein_name = data['sseqid']
        if protein_name not in protein_peptide_map:
            protein_peptide_map[protein_name] = []
        protein_peptide_map[protein_name].append(peptide)

    with open(fasta_file_path, 'r') as fasta_file:
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

    for peptide, data in peptide_map.items():
        if 'ssequence' not in data:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")

    return peptide_map

def output_peptide_data_for_reset(comet_map, casanovo_map, diamond_map, decoy_prefix):
    peptides = set(comet_map.keys()) | set(casanovo_map.keys())
    missing_peptides = peptides - set(diamond_map.keys())

    if missing_peptides:
        warning_message = f"Warning: The following {len(missing_peptides)} peptides are missing from diamond_map and will be skipped: " + ", ".join(missing_peptides)
        print(warning_message, file=sys.stderr)

    casanovo_charges = set(data['charge'] for data in casanovo_map.values())
    comet_charges = set(data['charge'] for data in comet_map.values())

    header = [
        "SpecId", "Label", "ScanNr", "database_peptide_length", "diamond_bitscore", "diamond_perc_identity",
        "casanovo_num_spectra", "casanovo_best_score", "casanovo_ppm_error", "casanovo_num_peptidoforms",
        "comet_num_spectra", "comet_n_tryptic", "comet_c_tryptic", "comet_best_score", "comet_ppm_error", "comet_num_peptidoforms",
        "combined_rank_score"
    ]
    # header.extend(f"casanovo_charge{charge}" for charge in casanovo_charges)
    # header.extend(f"comet_charge{charge}" for charge in comet_charges)
    header.extend(["Peptide", "Proteins"])
    print("\t".join(header))

    scan_nr = 1
    for peptide in peptides:
        if peptide not in diamond_map:
            continue

        diamond_data = diamond_map[peptide]
        casanovo_data = casanovo_map.get(peptide, {})
        comet_data = comet_map.get(peptide, {})

        spec_id = peptide + '_1'
        label = -1 if diamond_data.get('sseqid', '').startswith(decoy_prefix) or comet_data.get('is_decoy') == '1' else 1

        database_peptide_length = int(diamond_data.get('send', 0)) - int(diamond_data.get('sstart', 0)) + 1
        diamond_bitscore = diamond_data.get('bitscore', 0)
        diamond_perc_identity = diamond_data.get('pident', 0)

        casanovo_num_spectra = casanovo_data.get('num_spectra', 0)
        casanovo_best_score = casanovo_data.get('search_engine_score[1]', 0)
        casanovo_ppm_error = casanovo_data.get('mz_ppm_error', 0)
        casanovo_num_peptidoforms = casanovo_data.get('num_peptidoforms', 0)

        comet_num_spectra = comet_data.get('num_spectra', 0)
        comet_n_tryptic = comet_data.get('tryptic_n', 0)
        comet_c_tryptic = comet_data.get('tryptic_c', 0)
        comet_best_score = comet_data.get('e-value', 0)
        comet_ppm_error = comet_data.get('mz_ppm_error', 0)
        comet_num_peptidoforms = comet_data.get('num_peptidoforms', 0)

        combined_rank_score = str(4 - float(comet_data.get('rank_score', 2)) - float(casanovo_data.get('rank_score', 2)))

        ssequence = diamond_data.get('ssequence')
        if ssequence is None:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")

        proteins = diamond_data.get('sseqid', '')

        row = [
            spec_id, label, scan_nr, database_peptide_length, diamond_bitscore, diamond_perc_identity,
            casanovo_num_spectra, casanovo_best_score, casanovo_ppm_error, casanovo_num_peptidoforms,
            comet_num_spectra, comet_n_tryptic, comet_c_tryptic, comet_best_score, comet_ppm_error, comet_num_peptidoforms,
            combined_rank_score
        ]

        # for charge in casanovo_charges:
        #     row.append(1 if casanovo_data.get('charge') == charge else 0)

        # for charge in comet_charges:
        #     row.append(1 if comet_data.get('charge') == charge else 0)

        row.extend([ssequence, proteins])

        print("\t".join(str(value) for value in row))
        scan_nr += 1

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python build_reset_input.py <comet_results_file> <casanovo_results_file> <diamond_results_file> <fasta_file> <decoy_prefix>")
        sys.exit(1)

    comet_results_file = sys.argv[1]
    casanovo_results_file = sys.argv[2]
    diamond_results_file = sys.argv[3]
    fasta_file = sys.argv[4]
    decoy_prefix = sys.argv[5]

    comet_map = read_peptide_file(comet_results_file)
    casanovo_map = read_peptide_file(casanovo_results_file)
    diamond_map = read_diamond_file(diamond_results_file)
    diamond_map = augment_peptide_map(diamond_map, fasta_file)

    output_peptide_data_for_reset(comet_map, casanovo_map, diamond_map, decoy_prefix)