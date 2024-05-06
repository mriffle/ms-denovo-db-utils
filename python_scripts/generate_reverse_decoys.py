import sys
import argparse

def reverse_sequence(sequence):
    return sequence[::-1]

def print_entry(header, sequence, decoy_prefix):
    sequence = sequence.rstrip('*')
    print(f"{header}\n{sequence}")
    print(f">{decoy_prefix}{header[1:]}\n{reverse_sequence(sequence)}")

def process_fasta_file(input_file, decoy_prefix):
    header = ""
    sequence = ""
    with open(input_file, 'r') as file_in:
        for line in file_in:
            if line.startswith('>'):
                if header and sequence:
                    print_entry(header, sequence, decoy_prefix)
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
    if header and sequence:
        print_entry(header, sequence, decoy_prefix)

def main():
    parser = argparse.ArgumentParser(description='Process FASTA file and generate decoy sequences.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('--decoy_prefix', default='DECOY_', help='Decoy prefix (default: DECOY_)')

    args = parser.parse_args()
    input_file = args.input_file
    decoy_prefix = args.decoy_prefix

    process_fasta_file(input_file, decoy_prefix)

if __name__ == '__main__':
    main()
