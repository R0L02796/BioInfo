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

Para correr el programa

`pipenv run gb_to_fasta.py [-h] -i I [-o O]`

**Argumentos:**

- `-i I, --in-file I`: Archivo de entrada GeneBank con extensión .gb
- `-o O, --out-dir O`: Directorio de salida para las secuencias FASTA
    - default: `./out/fasta/`

Para este ejercicio ya se encuentra disponible el archivo GeneBank
correspondiente al gen HTT en `data/genebank/HTT_NM_001388492_iso1.gb`.

El programa toma el archivo en formato GeneBank y analiza las secuencias de
nucleótidos en sus registros. Luego, para cada una de ellas escribe un
archivo en formato FASTA como `<id>.fasta`.

También se encuentra un archivo `data/genebank/multi_record.gb` con 2 registros
para probar que el programa funciona con varios registros.


#### Ejercicio 2

Para correr el programa

`pipenv run blast.py [-h] -i I [-o O] -m {nucleotide,protein}`


**Argumentos:**

- `-i I, --in-dir I`: Directorio de entrada con secuencias FASTA
- `-o O, --out-dir O`: Directorio de salida para los reportes BLAST de cada secuencia
    - default: `./out/blast/`

El programa este se conecta con el API del NCBI y realiza el BLAST para cada
secuencia, escribiendo el resultado del informe en un archivo en formato texto.

Se puede correr el archivo utilizando como entrada el FASTA del  ejercicio
anterior `blast.py -i out/fasta -m nucleotide` (asumiendo que se corrió con la
salida default) para probar.

**Análisis de los resultados**

El BLAST encuentra la el gen NM_001388492.1 como la más parecida, con un
match del 100%. Esto es de esperarse, ya que esta es la variante del mRNA HTT que
levantó del archivo GenBank. Como segunda opción encuentra NM_002111.8 que es
otra variante, en este caso con 6 gaps y 0 mutaciones. Como tercera opción
encuentra la secuencia del gen de referencia que también es muy parecida, con
10 gaps nada más.

Varios de los resultados obtenidos ademas corresponden a primates lo cual es
razonable teniendo en cuenta que el genoma humano y el de los primates son muy
parecidos. Los siguientes resultados son de chimpancés y bonobos que no solo
entre sí como especie son muy similares, sino que son las especies vivas más
relacionadas a los humanos, lo cual explica también el score tan alto para
estos genes.
