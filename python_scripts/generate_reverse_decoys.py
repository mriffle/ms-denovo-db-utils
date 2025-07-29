import sys
import argparse
import gzip

def reverse_sequence(sequence):
    return sequence[::-1]

def print_entry(header, sequence, decoy_prefix, output_file):
    sequence = sequence.rstrip('*')
    output_file.write(f"{header}\n{sequence}\n")
    output_file.write(f">{decoy_prefix}{header[1:]}\n{reverse_sequence(sequence)}\n")

def is_gzipped(filename):
    """Check if file is gzipped by reading the magic number"""
    try:
        with open(filename, 'rb') as f:
            magic = f.read(2)
            return magic == b'\x1f\x8b'
    except IOError:
        return False

def process_fasta_file(input_file, decoy_prefix):
    # Determine if input is gzipped
    input_is_gzipped = is_gzipped(input_file)
    
    # Set up input file handle
    if input_is_gzipped:
        file_in = gzip.open(input_file, 'rt')
    else:
        file_in = open(input_file, 'r')
    
    # Set up output file handle
    if input_is_gzipped:
        output_file = gzip.open(sys.stdout.buffer, 'wt')
    else:
        output_file = sys.stdout
    
    try:
        header = ""
        sequence = ""
        
        for line in file_in:
            if line.startswith('>'):
                if header and sequence:
                    print_entry(header, sequence, decoy_prefix, output_file)
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        
        # Process the last entry
        if header and sequence:
            print_entry(header, sequence, decoy_prefix, output_file)
            
    finally:
        file_in.close()
        if input_is_gzipped:
            output_file.close()

def main():
    parser = argparse.ArgumentParser(description='Process FASTA file and generate decoy sequences.')
    parser.add_argument('input_file', help='Input FASTA file (can be gzipped)')
    parser.add_argument('--decoy_prefix', default='DECOY_', help='Decoy prefix (default: DECOY_)')

    args = parser.parse_args()
    input_file = args.input_file
    decoy_prefix = args.decoy_prefix

    process_fasta_file(input_file, decoy_prefix)

if __name__ == '__main__':
    main()