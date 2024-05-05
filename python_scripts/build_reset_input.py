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

def read_glsearch36_file(file_path):
    column_headers = [
        'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
        'sstart', 'send', 'evalue', 'bitscore'
    ]
    expected_columns = len(column_headers) + 1
    peptide_map = {}

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            columns = line.strip().split('\t')

            if len(columns) != expected_columns:
                raise ValueError(f"Line {line_num}: Expected {expected_columns} columns, but found {len(columns)}")

            peptide_sequence = columns[0]
            peptide_data = {column_headers[i]: columns[i+1] for i in range(len(column_headers))}
            peptide_map[peptide_sequence] = peptide_data

    return peptide_map

def parse_glsearch36_directory(directory_path):
    peptide_map = {}

    for file_name in os.listdir(directory_path):
        if file_name.endswith("gl.txt"):
            file_path = os.path.join(directory_path, file_name)
            file_peptide_map = read_glsearch36_file(file_path)

            for peptide_sequence in file_peptide_map:
                if peptide_sequence in peptide_map:
                    raise ValueError(f"Duplicate peptide sequence found: {peptide_sequence}")

            peptide_map.update(file_peptide_map)

    return peptide_map

def augment_peptide_map(peptide_map, fasta_file_path):
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

def output_peptide_data_for_reset(comet_map, casanovo_map, glsearch_map, decoy_prefix):
    peptides = set(comet_map.keys()) | set(casanovo_map.keys())
    missing_peptides = peptides - set(glsearch_map.keys())

    if missing_peptides:
        warning_message = "Warning: The following peptides are missing from glsearch_map and will be skipped: " + ", ".join(missing_peptides)
        print(warning_message, file=sys.stderr)

    casanovo_charges = set(data['charge'] for data in casanovo_map.values())
    comet_charges = set(data['charge'] for data in comet_map.values())

    header = [
        "SpecId", "Label", "ScanNr", "database_peptide_length", "glsearch_evalue", "glsearch_perc_identity",
        "casanovo_num_spectra", "casanovo_best_score", "casanovo_ppm_error",
        "comet_num_spectra", "comet_n_tryptic", "comet_c_tryptic", "comet_best_score", "comet_ppm_error"
    ]
    header.extend(f"casanovo_charge{charge}" for charge in casanovo_charges)
    header.extend(f"comet_charge{charge}" for charge in comet_charges)
    header.extend(["Peptide", "Proteins"])
    print("\t".join(header))

    scan_nr = 1
    for peptide in peptides:
        if peptide not in glsearch_map:
            continue

        glsearch_data = glsearch_map[peptide]
        casanovo_data = casanovo_map.get(peptide, {})
        comet_data = comet_map.get(peptide, {})

        spec_id = peptide
        label = -1 if glsearch_data.get('sseqid', '').startswith(decoy_prefix) else 1

        database_peptide_length = int(glsearch_data.get('send', 0)) - int(glsearch_data.get('sstart', 0)) + 1
        glsearch_evalue = glsearch_data.get('evalue', 0)
        glsearch_perc_identity = glsearch_data.get('pident', 0)

        casanovo_num_spectra = casanovo_data.get('num_spectra', 0)
        casanovo_best_score = casanovo_data.get('search_engine_score[1]', 0)
        casanovo_ppm_error = casanovo_data.get('mz_ppm_error', 0)

        comet_num_spectra = comet_data.get('num_spectra', 0)
        comet_n_tryptic = comet_data.get('tryptic_n', 0)
        comet_c_tryptic = comet_data.get('tryptic_c', 0)
        comet_best_score = comet_data.get('e-value', 0)
        comet_ppm_error = comet_data.get('mz_ppm_error', 0)

        ssequence = glsearch_data.get('ssequence')
        if ssequence is None:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")

        proteins = glsearch_data.get('sseqid', '')

        row = [
            spec_id, label, scan_nr, database_peptide_length, glsearch_evalue, glsearch_perc_identity,
            casanovo_num_spectra, casanovo_best_score, casanovo_ppm_error,
            comet_num_spectra, comet_n_tryptic, comet_c_tryptic, comet_best_score, comet_ppm_error
        ]

        for charge in casanovo_charges:
            row.append(1 if casanovo_data.get('charge') == charge else 0)

        for charge in comet_charges:
            row.append(1 if comet_data.get('charge') == charge else 0)

        row.extend([ssequence, proteins])

        print("\t".join(str(value) for value in row))
        scan_nr += 1

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python build_reset_input.py <comet_results_file> <casanovo_results_file> <glsearch36_directory> <fasta_file> <decoy_prefix>")
        sys.exit(1)

    comet_results_file = sys.argv[1]
    casanovo_results_file = sys.argv[2]
    glsearch36_directory = sys.argv[3]
    fasta_file = sys.argv[4]
    decoy_prefix = sys.argv[5]

    comet_map = read_peptide_file(comet_results_file)
    casanovo_map = read_peptide_file(casanovo_results_file)
    glsearch_map = parse_glsearch36_directory(glsearch36_directory)
    glsearch_map = augment_peptide_map(glsearch_map, fasta_file)

    output_peptide_data_for_reset(comet_map, casanovo_map, glsearch_map, decoy_prefix)
