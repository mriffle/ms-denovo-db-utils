import sys

def reverse_sequence(sequence):
    return sequence[::-1]

def process_fasta_file(input_file):
    header = ""
    sequence = ""

    with open(input_file, 'r') as file_in:
        for line in file_in:
            if line.startswith('>'):
                if header and sequence:
                    sequence = sequence.rstrip('*')
                    print(f"{header}\n{sequence}")
                    print(f">DECOY_{header[1:]}\n{reverse_sequence(sequence)}")
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()

        if header and sequence:
            sequence = sequence.rstrip('*')
            print(f"{header}\n{sequence}")
            print(f">DECOY_{header[1:]}\n{reverse_sequence(sequence)}")

# Check if the input file is provided as a command-line argument
if len(sys.argv) < 2:
    print("Please provide the input FASTA file as a command-line argument.")
    sys.exit(1)

input_file = sys.argv[1]
process_fasta_file(input_file)
