# file stream script to analyze CpA methylation
from asyncore import read
import csv
import os
import glob
import argparse


def txtToBedFile(file, args):
    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        # design output file name based off the input file (subsitute .txt for .bed)
        # output_name = f"{tab_file.name.split('.')[0]}.bed"
        output_name = tab_file.name.replace('.txt', '.bed')

        # if the output file does not already exist proceed
        if not os.path.exists(output_name):
            writeBedFile(output_name, reader, args)


def writeBedFile(output_name, reader, args):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chr', 'start', 'end', 'strand']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        # write header to output file if the header flag is true
        if args.header:
            writer.writeheader()

        # loop through the lines read in
        for line in reader:
            # a list of columns from the input file to be removed
            keys_to_remove = ['chrBase', 'coverage', 'freqC', 'freqT']

            # remove undesired columns from read in line
            for k in keys_to_remove:
                line.pop(k, None)

            # rename column 'base' to 'start'
            line['start'] = line.pop('base')

            # generate the end column by adding 1 to the start column value
            line['end'] = int(line.get('start')) + 1

            # covert the 'F' and 'R' strand notation for '+' and '-'
            if line.get('strand') == 'F':
                line['strand'] = '+'
            else:
                line['strand'] = '-'

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
    parser.add_argument('-H', '--header', action='store_true', help='if flag is provided a header line will be written to the output file')
    args = parser.parse_args()

    # navigate to the specified directory
    os.chdir(args.dir)

    # glob all files in the current directory following the specified pattern: example '*_CHH.txt'
    for file in glob.glob(args.pattern):
        txtToBedFile(file, args)


if __name__ == "__main__":
    main()
