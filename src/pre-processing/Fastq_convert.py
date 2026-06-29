import gzip
import argparse

from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Convert a FASTQ file to FASTA format.")
    parser.add_argument("input_file", help="The name of the input FASTQ file (gzip compressed).")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = f"{input_file[:-6]}.fasta"

    with gzip.open(input_file, "rt") as fastq_handle:
        with open(output_file, "w") as fasta_handle:
            SeqIO.convert(fastq_handle, "fastq", fasta_handle, "fasta")

    print(f"[-] Conversion complete: saved as {output_file}.")

if __name__ == "__main__":
    main()