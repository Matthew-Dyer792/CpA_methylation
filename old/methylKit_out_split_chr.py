# file stream script to analyze CpA methylation
import csv
from pathlib import Path
import time


FILE_PATH = '/research/project/shared/benoukraf_lab/matthew/DNMT/WGBS_on_sorted_PSCs-6_samples/meth_raw_CHH/DNMT21-1-T22-G1-1_CHH.txt'


def read_methylkit_out(file, out_dir):
    # design output file name based off the input file (substitute .txt for .base.CpA.txt)
    original_name = Path(file).name

    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        for line in reader:
            # call the 'write_bed' function
            write_file(out_dir, original_name, line)


def write_file(out_dir, original_name, line):
    # get current lines chromosome
    chr = line.get('chr')
    # add chromosome name to the output file
    output_name = original_name.replace('.txt', f".{chr}.txt")
    # merge the out directory with the output name
    output_name = out_dir / output_name
    
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chrBase', 'chr', 'base', 'strand', 'coverage', 'freqC', 'freqT']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        # write modified line to the output file in append
        writer.writerow(line)


def main():
    file = Path(FILE_PATH)
    out_dir = Path(FILE_PATH.replace('.txt', '_chr'))

    # creates path and parents if not already existing
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    start = time.time()

    read_methylkit_out(file, out_dir)

    end = time.time()
    print(end - start)


if __name__ == "__main__":
    main()
