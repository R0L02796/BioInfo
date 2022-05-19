# Introducción a la Bioinformática
<br />
<p align="center">
  <a href="https://github.com/R0L02796/BioInfo">
    <img src="assets/logo.png" alt="Logo" width="234" height="115">
  </a>
</p>

<h3 align="center">Enfermedad de Huntington</h3>

## Tabla de Contenidos

- [Introducción a la Bioinformática - Trabajo Práctico](#introducción-a-la-bioinformática)
  - [Tabla de Contenidos](#tabla-de-contenidos)
  - [Introducción](#introduccion)
    - [Requisitos](#requisitos)
    - [Instalación](#instalación)
    - [Ejecución](#ejecución)
      -[Ejercicio 1](#ejercicio-1)

## Introducción

Trabajo práctico para la materia Introducción a la Bioinformática con el
objetivo de adquirir habilidades básicas en el campo de la Bioinformática. Para
este trabajo práctico se decidió investigar el gen HTT, que está relacionado
con la enfermedad de Huntington.

### Requisitos

- Python3
- pip3
- [pipenv](https://pypi.org/project/pipenv/)

### Instalación

Parado en la carpeta raíz del proyecto hacer

```sh
pipenv install
```

Para instalar las dependencias necesarias en el ambiente virtual.

### Ejecución

#### Ejercicio 1

`python3 gb_to_fasta.py [-h] -i I [-o O]`

**Argumentos:**

- `-i I, --in_file I`: Archivo de entrada GeneBank con extensión .gb
- `-o O, --out_dir O`: Directorio de salida para las secuencias FASTA
    - default: `./out/fasta/`

Para este ejercicio ya se encuentra disponible el archivo GeneBank
correspondiente al gen HTT en `data/genebank/HTT_NM_001388492_iso1.gb`.

El programa toma el archivo en formato GeneBank y analiza las secuencias de
nucleótidos en sus registros. Luego, para cada una de ellas escribe un
archivo en formato FASTA como `<id>.fasta`.



