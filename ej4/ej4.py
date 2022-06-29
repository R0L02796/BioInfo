import argparse
import os
from Bio import Entrez
from Bio.Blast import NCBIXML

DEFAULT_BLAST_DIR = './out/blast-output/'
FASTA_LINE_LENGTH = 70
FASTA = "fasta"

parser = argparse.ArgumentParser(description='Runs BLAST output')
parser.add_argument('-i', '--in-file', help='', metavar="I", type=str, required=True)
parser.add_argument('-o', '--out-file', help='Pattern analisis outfile', metavar="O",  type=str, required=False, default=DEFAULT_BLAST_DIR+"ej4.txt")
parser.add_argument('-m', '--mode', help='Choose between nucleotide or protein BLAST', choices=('nucleotide', 'protein'), type=str, required=True)
parser.add_argument('-p', '--pattern', help='Pattern', type=str, required=True)
args = parser.parse_args()

blast_filename = args.in_file
faa_filename = args.out_file
mode = args.mode

# Groundwork for more error handling
ERROR_TYPE = None
ERROR_MSGS = {
    "IN_FILE_NOT_VALID": "El archivo de entrada no es de tipo blast output (xml)",
    "OUT_FILE_NOT_VALID": "El archivo de salida no es de tipo txt",
}

try:
    if(args.in_file.split(".")[-1] != "out"):
        ERROR_TYPE = "IN_FILE_NOT_VALID"
        raise Exception
    if(args.out_file.split(".")[-1] != "txt"):
        ERROR_TYPE = "OUT_FILE_NOT_VALID"
        raise Exception

except Exception:
    print(ERROR_MSGS[ERROR_TYPE])
    exit(1)

out_dir = os.path.dirname(args.out_file)
if(out_dir != ''):
    # Generate directory if necessary for output file
    os.makedirs(out_dir, exist_ok=True)


in_file = open(blast_filename)
save_file = open(faa_filename, "w")
record = NCBIXML.parse(in_file)

Entrez.email = "@"
for rec in record:
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if(args.pattern.lower() in alignment.title.lower()):
                save_file.write('****Alignment****\n')
                save_file.write('-sequence: ' + alignment.title + '\n')
                save_file.write('-accession: ' + alignment.accession + '\n')
                save_file.write('-length: ' + str(alignment.length) + '\n')
                save_file.write('-score: ' + str(hsp.score) + '\n')
                save_file.write('-gaps: ' + str(hsp.gaps) + '\n')
                if type == "protein":
                    save_file.write(Entrez.efetch(
                        db="protein", rettype="fasta", id=alignment.accession).read())
                else:
                    save_file.write(Entrez.efetch(
                        db="nuccore", rettype="fasta", id=alignment.accession).read())

save_file.close()
