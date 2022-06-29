import docker
import sys
import os

client = docker.from_env()

if len(sys.argv) != 3:
    print("Wrong arguments amount")
    exit(1)

faa_filename = sys.argv[1]
out_filename = sys.argv[2]

if faa_filename.split(".")[1] != "fasta" or out_filename.split(".")[1] != "patmatmotifs":
    print("Wrong format for input")
    exit(1)

client.containers.run("pegi3s/emboss:latest", "transeq -sequence " + faa_filename + " -outseq EJ5/ej5.transeq.fa")
client.containers.run("pegi3s/emboss:latest", "patmatmotifs -sequence EJ5/ej5.transeq.fa -outfile " + out_filename)
