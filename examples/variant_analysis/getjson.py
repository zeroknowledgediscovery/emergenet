from Bio import SeqIO
import json
import sys

def extract_sequence(fasta_file):
    """Extracts and concatenates sequences from a FASTA file."""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return ''.join(sequences)

def main(file_prefix):
    # Constructing file paths
    ha_fasta = f"{file_prefix}_ha.fasta"
    na_fasta = f"{file_prefix}_na.fasta"

    # Extract sequences
    ha_sequence = extract_sequence(ha_fasta)
    na_sequence = extract_sequence(na_fasta)

    # Output in JSON format
    output = {"id":file_prefix.replace(':','/'),
              "HA": ha_sequence, "NA": na_sequence}
    print(json.dumps(output, indent=4))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py file_prefix")
        sys.exit(1)
    
    file_prefix = sys.argv[1]
    main(file_prefix)
