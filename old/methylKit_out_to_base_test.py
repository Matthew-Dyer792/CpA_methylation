# file stream script to analyze CpA methylation
import csv
from pathlib import Path
import time
import pysam


# conversion rules for reverse strand
REVERSE_ALPHABET = {
    'a': 't', 'A': 'T',
    't': 'a', 'T': 'A',
    'c': 'g', 'C': 'G',
    'g': 'c', 'G': 'c'
}
# list of chromosomes
CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8'
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


def read_methylkit_out(file, referenceGenome):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    output_name = Path(Path(file).name.replace('.txt', '.CpA.bed'))

    # if the output file exists delete it and proceed
    try:
        output_name.unlink(missing_ok=True)
    except OSError:
        pass

    chh_set = set([])

    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        for line in reader:
            # call the 'write_bed' function
            write_file(output_name, line, referenceGenome, chh_set)

    return chh_set


def write_file(output_name, line, referenceGenome, chh_set):
    # define the chromosome of the current line
    chr = line.get('chr')

    if chr in CHROMOSOMES:
        # define the genoimic location for the nucleotide of the current line
        base = int(line.get('base'))

        # open output file in append mode
        with open(output_name, 'a') as output_file:
            # set the names of the fields to be written to the new file
            field_names = ['chr', 'start', 'end', 'strand', 'coverage', 'freqC', 'freqT', 'nucleotides']

            # create line writer
            writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')
            
            if line.get('strand') == 'R':
                # get the three bases after the C nucleotide
                nucleotide = referenceGenome.fetch(chr, base - 3, base)

                chh_set.add(nucleotide.upper())

                # reverse letters to match forward direction
                nucleotide = nucleotide[::-1]
                # transform the letters to their forward compliment
                letters = [REVERSE_ALPHABET.get(letter) for letter in list(nucleotide)]

                if letters is not None:
                    # join the letters back into a str
                    nucleotide = "".join(letters)
                else:
                    chh_set.add(f"error on this set: {nucleotide.upper()}")
            else:
                # get the three bases after the C nucleotide
                nucleotide = referenceGenome.fetch(chr, base - 1, base + 2)

                chh_set.add(nucleotide.upper())

            # covert the 'F' and 'R' strand notation for '+' and '-'
            if line.get('strand') == 'F':
                line['strand'] = '+'
            else:
                line['strand'] = '-'

            # only write out sites with a CpA context 
            if nucleotide.upper().startswith('CA')
                # write modified line to the output file in append
                writer.writerow(line)


def main():
    # open reference genome fasta file location
    referenceGenome = pysam.Fastafile(Path('/research/project/shared/benoukraf_lab/matthew/data/GRCh38_UCSC/hg38.fa'))

    file = Path('/research/project/shared/benoukraf_lab/matthew/DNMT/WGBS_on_sorted_PSCs-6_samples/meth_raw_CHH/DNMT21-1-T22-G1-1_CHH.txt')

    start = time.time()

    chh_set = read_methylkit_out(file, referenceGenome)

    end = time.time()
    print(end - start)
    print(chh_set)

    referenceGenome.close()



if __name__ == "__main__":
    main()
