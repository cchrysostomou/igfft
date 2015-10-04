# igfft

FFT implementation for NGS antibody sequence alignments

## Description

This program is used for aligning NGS data to germline antibody genes. It utilizes the FFTW program to perform local 
alignments of a read to a provided list of possible genes. On a single 2.7 GHz processor, it can align 1 million sequences to germline genes
in 20-30 minutes.

## General usage

Inputs:
The minimum requirement for using the program is to define an input file of NGS data (FASTA,FASTQ,TAB format) and to define the location of at least
one database file (V or J germlines). 

Output:
The output of the program is a single tab delimited file summarizing (1) V and/or J hits to read, (2) full length antibody region, (3) annotated antibody regions (FR1,FR2,FR3,CDR1,CDR2,CDR3-FR4), (4) Alignments to the top V and or J germline gene

## Examples
Running program from a terminal

### In this example we run the program by defining the input FASTQ file, the location of the V and J germline files, and an output file 
	binary/igfft sample/Demo_2.query.fastq -i FASTQ -v germlines/HomosapiensIGH_IGK_IGL_V.txt -j HomosapiensIGH_IGK_IGL_J.txt

### In this next example we will also define the number of top v hits we want returned (2 in this example)
	binary/igfft sample/Demo_2.query.fastq -i FASTQ -v germlines/HomosapiensIGH_IGK_IGL_V.txt -j HomosapiensIGH_IGK_IGL_J.txt -o sample/resultfile.txt -vnum_hits 2

### All possible default parameters in the program can be found using the following command
binary/igfft --defaults

## Dependencies
Installation of the FFTW library package is required to run this program. 
The binary provided was compiled on a linux machine. 

## Compiling source code
The following command can be used to compile source code:
g++ MainProgram.cpp DnaFunctions.cpp FFTFunctions.cpp ReadFiles.cpp GermlineCluster.cpp QueryAlignment.cpp SWAlignment.cpp FFTAlign.cpp -std=c++11 -lfftw3 -O



