#!/usr/bin/env bash
# Usage: split_fasta.sh <fasta_file> <number_of_parts>

# Check if sufficient arguments were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <number_of_parts>"
    exit 1
fi

fasta_file=$1  # First parameter: FASTA file name
user_defined_parts=$2  # Second parameter: Number of parts

# Validate the input parameters
if [ ! -f "$fasta_file" ]; then
    echo "Error: File '$fasta_file' does not exist."
    exit 1
fi

if ! [[ $user_defined_parts =~ ^[0-9]+$ ]]; then
    echo "Error: Number of parts must be an integer."
    exit 1
fi

# Count the number of sequences in the file.
total_sequences=$(grep -c "^>" "$fasta_file")

# Determine the minimum of the total sequences and the user-defined number of parts.
N=$((total_sequences < user_defined_parts ? total_sequences : user_defined_parts))

# Split the FASTA file while keeping sequences intact.
awk -v N=$N 'BEGIN {file_number = 1; sequence_number = 0; file_name = sprintf("query_part%d.fasta", file_number);}
  /^>/ {
    if (file_number < N && sequence_number >= int((NR-1)/N) + 1) {
      close(file_name);
      file_number++;
      file_name = sprintf("query_part%d.fasta", file_number);
      sequence_number = 0;
    }
  }
  { print >> file_name; if (/^>/) sequence_number++; }
' "$fasta_file"
