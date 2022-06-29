import docker
import shutil
import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--in-file', help='Input FASTA file', metavar="I", type=str, required=True)
parser.add_argument('-o', '--out-file', help='Output file', metavar="O",  type=str, required=False, default="ej5.patmatmotifs")
args = parser.parse_args()

# Groundwork for more error handling
ERROR_TYPE = None
ERROR_MSGS = {
    "IN_FILE_NOT_VALID": "El archivo de entrada no es de tipo fasta",
    "OUT_FILE_NOT_VALID": "El archivo de salida no es de tipo matmatmotifs",
}

try:
    if(args.in_file.split(".")[1] != "fasta"):
        ERROR_TYPE = "IN_FILE_NOT_VALID"
        raise Exception
    if(args.out_file.split(".")[1] != "patmatmotifs"):
        ERROR_TYPE = "OUT_FILE_NOT_VALID"
        raise Exception

except Exception:
    print(ERROR_MSGS[ERROR_TYPE])
    exit(1)

out_dir = os.path.dirname(args.out_file)
if(out_dir != ''):
    # Generate directory if necessary for output file
    os.makedirs(out_dir, exist_ok=True)

shutil.copyfile(args.in_file, "data/file.fasta")
client = docker.from_env()
print("Running EMBOSS operations")
cmd = 'sh -c "transeq -sequence /data/file.fasta -outseq /data/tmp.fa && prosextract -prositedir /data/prosite && patmatmotifs -sequence /data/tmp.fa -outfile /data/out.patmatmotifs"'
client.containers.run("pegi3s/emboss:latest",cmd, auto_remove=True, volumes=[f"{os.getcwd()}/data:/data"])
shutil.copyfile("data/out.patmatmotifs", args.out_file)


