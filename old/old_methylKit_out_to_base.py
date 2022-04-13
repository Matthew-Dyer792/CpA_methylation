# file stream script to analyze CpA methylation
import csv
import os
from pathlib import Path
import glob
import argparse
import pysam


# # open genome fasta file location
# GENOME = pysam.Fastafile('/home/matthew/data/hg38.fa')
# conversion rules for reverse strand
REVERSE_ALPHABET = {
    'a': 't', 'A': 'T',
    't': 'a', 'T': 'A',
    'c': 'g', 'C': 'G',
    'g': 'c', 'G': 'c'
}


def read_methylkit_out(file, args):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    output_name = args.outDir + '/' + Path(file).name.replace('.txt', '.base.CpA.txt')

    # if the output file exists delete it and proceed
    try:
        os.remove(output_name)
    except OSError:
        pass

    # write header to output file if the header flag is true
    if args.header:
        write_header(output_name)

    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        for line in reader:
            # call the 'write_bed' function
            write_file(output_name, line, args)


def write_header(output_name):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chrBase', 'chr', 'base', 'strand', 'coverage', 'freqC', 'freqT', 'nucleotide']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        writer.writeheader()


def write_file(output_name, line, args):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chrBase', 'chr', 'base', 'strand', 'coverage', 'freqC', 'freqT', 'nucleotide']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        chr = line.get('chr')
        base = int(line.get('base'))
        
        if line.get('strand') == 'R':
            # get the three bases after the C nucleotide
            nucleotide = args.referenceGenome.fetch(chr, base - 3, base)
            # reverse letters to match forward direction
            nucleotide = nucleotide[::-1]
            # transform the letters to their forward compliment
            letters = [REVERSE_ALPHABET.get(letter) for letter in list(nucleotide)]
            # join the letters back into a str
            nucleotide = "".join(letters)
        else:
            # get the three bases after the C nucleotide
            nucleotide = args.referenceGenome.fetch(chr, base - 1, base + 2)

        line['nucleotide'] = nucleotide

        # only write CpA bases
        if nucleotide.lower().startswith('ca'):
            # write modified line to the output file in append
            writer.writerow(line)


def main():
    # setup command line arguments to be passed to the script
    parser = argparse.ArgumentParser(
        prog='MethylKit out to BED file',
        description='A program to convert the output of the methylKit CHH analysis to a BED file to be passed to bedtools function "getFasta"'
    )
    parser.add_argument('-d', '--dir', nargs='?', required=True, help='directory containing methylKit\'s output')
    parser.add_argument('-p', '--pattern', nargs='?', required=True, help='file pattern to glob output files')
    parser.add_argument('-o', '--outDir', nargs='?', help='directory for the output file(s)')
    parser.add_argument('-r', '--referenceGenome', nargs='?',required=True,  help='path to your reference genome fasta file')
    parser.add_argument('-H', '--header', action='store_true', help='if flag is provided a header line will be written to the output file')
    args = parser.parse_args()

    # directory containing the raw CHH files
    args.dir = os.path.abspath(args.dir)

    if args.outDir is not None:
        args.outDir = os.path.abspath(args.outDir)
    else:
        # default output if args.outDir not given
        args.outDir = args.dir + '/cpa_context'

    # creates path and parents if not already existing
    Path(args.outDir).mkdir(parents=True, exist_ok=True)

    # open reference genome fasta file location
    args.referenceGenome = pysam.Fastafile(os.path.abspath(args.referenceGenome))

    # glob all files in the current directory following the specified pattern: example '*_CHH.txt'
    for file in glob.glob(f"{args.dir}/{args.pattern}"):
        read_methylkit_out(file, args)

    args.referenceGenome.close()


if __name__ == "__main__":
    main()
