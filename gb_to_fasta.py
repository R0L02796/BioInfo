import argparse
import os
from Bio import SeqIO

DEFAULT_FASTA_FILEPATH = './out/fasta/'
FASTA_LINE_LENGTH = 70
GENBANK = "gb"

parser = argparse.ArgumentParser(description='Get FASTA sequence from Genebank file')
parser.add_argument('-i', '--in-file', help='Genebank input file', metavar="I", type=str, required=True)
parser.add_argument('-o', '--out-dir', help='FASTA output dir', metavar="O",  type=str, required=False, default=DEFAULT_FASTA_FILEPATH)
args = parser.parse_args()

gb_path = args.in_file
fasta_dir = args.out_dir

# Groundwork for more error handling

ERROR_TYPE = None
ERROR_MSGS = {
    "GB_EXT": "Archivo de entrada debería tener la extensión .gb",
}

try:
    if( gb_path.split(".")[1] != "gb"):
        ERROR_TYPE = "GB_EXT"
        raise Exception

except IndexError:
    print("El archivo de entrada debe tener una extensión")
    exit(1)
except Exception:
    print(ERROR_MSGS[ERROR_TYPE])
    exit(1)


# Checks if file is in root directory
if(fasta_dir != ''):
    # Generate directory if necessary for output file
    os.makedirs(fasta_dir, exist_ok=True)

for gb_content in SeqIO.parse(gb_path, GENBANK):
    # extracting header
    header = ">{id} {description}\n".format(id=gb_content.id, description=gb_content.description)
    # extracting sequence
    sequence = str(gb_content.seq)
    # making a list of strings with fasta fixed line length
    sequence_lines = [sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)]

    try:
        with open (f"{fasta_dir}{gb_content.id}.fasta", "w+") as fasta_file:
            # writing header to file 
            fasta_file.write(header)
            # writing each sequence line 
            fasta_file.writelines(["%s\n" % seq_line  for seq_line in sequence_lines])
            # adding \n at end of file
            fasta_file.write('\n')
    except:
        print("No se pudo escribir al directorio")
        exit(1)
