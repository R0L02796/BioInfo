from Bio import SeqIO

gb_path = "./sequence.gb"
fasta_path = "./sequence.fasta"

FASTA_LINE_LENGTH = 70
GENBANK = "gb" 


for gb_content in SeqIO.parse(gb_path, GENBANK):
    # extracting header
    header = ">{id} {description}\n".format(id=gb_content.id, description=gb_content.description)
    # extracting sequence
    sequence = str(gb_content.seq)
    # making a list of strings with fasta fixed line length
    sequence_lines = [sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)]
    # open fasta file
    fasta_file = open(fasta_path, "w+")
    # writing header to file 
    fasta_file.write(header)
    # writing each sequence line 
    fasta_file.writelines(["%s\n" % seq_line  for seq_line in sequence_lines])
    # adding \n at end of file
    fasta_file.write('\n')


fasta_file.close()