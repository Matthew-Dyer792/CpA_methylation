# file stream script to analyze CpA methylation
import csv
from pathlib import Path
import argparse
import pysam


# conversion rules for reverse strand
REVERSE_ALPHABET = {
    'a': 't', 'A': 'T',
    't': 'a', 'T': 'A',
    'c': 'g', 'C': 'G',
    'g': 'c', 'G': 'c'
}


def read_methylkit_out(args):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    output_name = args.outDir / args.file.name.replace('.txt', '.CpA.bed')

    write_header(output_name)

    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(args.file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        for line in reader:
            # call the 'write_bed' function
            write_file(output_name, line, args)


def write_header(output_name):
    # open output file in append mode
    with open(output_name, 'w') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chr', 'start', 'end', 'strand', 'coverage', 'freqC', 'freqT', 'nucleotide']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        writer.writeheader()


def write_file(output_name, line, args):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chr', 'start', 'end', 'strand', 'coverage', 'freqC', 'freqT', 'nucleotide']

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

            # assign start and end
            line['start'] = base - 2
            line['end'] = base + 1
        else:
            # get the three bases after the C nucleotide
            nucleotide = args.referenceGenome.fetch(chr, base - 1, base + 2)

            # assign start and end
            line['start'] = base
            line['end'] = base + 3

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
    parser.add_argument('-f', '--file', nargs='?', required=True, help='methylKit\'s output file')
    parser.add_argument('-o', '--outDir', nargs='?', help='directory for the output file(s)')
    parser.add_argument('-r', '--referenceGenome', nargs='?',required=True,  help='path to your reference genome fasta file')
    args = parser.parse_args()

    # directory to store the results
    args.outDir = Path(args.file.replace('.txt', '_chr'))

    # directory containing the raw CHH files
    args.file = Path(args.file)

    # creates path and parents if not already existing
    Path(args.outDir).mkdir(parents=True, exist_ok=True)

    # open reference genome fasta file location
    args.referenceGenome = pysam.Fastafile(Path(args.referenceGenome))

    read_methylkit_out(args)

    # close reference genome
    args.referenceGenome.close()


if __name__ == "__main__":
    main()
