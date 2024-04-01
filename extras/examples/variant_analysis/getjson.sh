#!/bin/bash

# Check if an argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 file_prefix"
    exit 1
fi

# File prefix from the argument
file_prefix=$1

# Constructing file paths
ha_fasta="${file_prefix}_ha.fasta"
na_fasta="${file_prefix}_na.fasta"

# Function to extract and concatenate sequences from a FASTA file
extract_sequence() {
    fasta_file=$1
    # Skip the first line (sequence ID), concatenate the remaining lines
    tail -n +2 "$fasta_file" | tr -d '\n'
}

# Extract sequences
ha_sequence=$(extract_sequence "$ha_fasta")
na_sequence=$(extract_sequence "$na_fasta")

# Output in JSON format
echo "{\"HA\": \"$ha_sequence\", \"NA\": \"$na_sequence\"}"
