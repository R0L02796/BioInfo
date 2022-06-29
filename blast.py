#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO
from Bio.Blast.NCBIWWW import qblast

DEFAULT_BLAST_DIR = './out/blast/'
FASTA_LINE_LENGTH = 70
FASTA = "fasta"

parser = argparse.ArgumentParser(description='Runs BLAST of FASTA sequences')
parser.add_argument('-i', '--in-dir', help='Directory with FASTA sequencies', metavar="I", type=str, required=True)
parser.add_argument('-o', '--out-dir', help='Directory for BLAST reports', metavar="O",  type=str, required=False, default=DEFAULT_BLAST_DIR)
parser.add_argument('-m', '--mode', help='Choose between nucleotide or protein BLAST', choices=('nucleotide', 'protein', 'nucleotide-xml', 'protein-xml'), type=str, required=True)
args = parser.parse_args()

fasta_dir = args.in_dir
blast_dir = args.out_dir
mode = args.mode

# Groundwork for more error handling
ERROR_TYPE = None
ERROR_MSGS = {
    "IN_DIR_NOT_VALID": "Directorio de entrada no es accesible",
    "NO_INPUT_FILES": "No se encuentran archivos fasta en el directorio"
}

try:
    if(os.path.isdir(fasta_dir) is False):
        ERROR_TYPE = "IN_DIR_NOT_VALID"
        raise Exception

    if(len([f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]) == 0):
        ERROR_TYPE = "NO_INPUT_FILES"
        raise Exception

except Exception:
    print(ERROR_MSGS[ERROR_TYPE])
    exit(1)

if(blast_dir != ''):
    # Generate directory if necessary for output file
    os.makedirs(blast_dir, exist_ok=True)


def run_blast(file):
    BLAST_MODE = {
        'nucleotide': lambda : qblast("blastn", "nt", record.format("fasta"), expect=0.05, format_type="Text"),
        'protein': lambda : qblast("blastn", "nr", record.format("fasta"), expect=0.05, format_type="Text"),
        'nucleotide-xml': lambda : qblast("blastn", "nt", record.format("fasta"), expect=0.05),
        'protein-xml': lambda : qblast("blastn", "nr", record.format("fasta"), expect=0.05),
    }

    for record in SeqIO.parse(file, FASTA):
        print(f"Running BLAST using NCBI API for sequence {record.id}...")
        result = BLAST_MODE[mode]()
        fp = os.path.join(blast_dir, f"{record.id}.out")
        print(f"Writing results to file: {fp}")
        with open(fp, "w+") as f:
            f.write(result.read())



for file in os.listdir(fasta_dir):
    if file.endswith(".fasta"):
        run_blast(os.path.join(fasta_dir, file))
    else:
        continue

