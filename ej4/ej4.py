import sys
from Bio import Entrez
from Bio.Blast import NCBIXML


if len(sys.argv) != 5:
    print("4 arguments expected")
    exit(1)

blast_filename = sys.argv[1]
faa_filename = sys.argv[2]

if blast_filename.split(".")[1] != "out":
    print("Input file must be in Blast format")
    exit(1)

pattern = sys.argv[3]
type = sys.argv[4]

if   type != 'nucleotide' and type != 'protein':
    print("Types not supported, try nucleotide or protein")
    exit(1)

in_file = open(blast_filename)
save_file = open(faa_filename, "w")

record = NCBIXML.parse(in_file)

Entrez.email = "@"

for rec in record:
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if(pattern.lower() in alignment.title.lower()):
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