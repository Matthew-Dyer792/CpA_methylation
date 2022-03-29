# file stream script to analyze CpA methylation
from asyncore import read
from copy import deepcopy
import csv
import os
import glob
import argparse


# file_line = 0


def txt_to_bed_file(file, args):
    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        # design output file name based off the input file (subsitute .txt for .bed)
        output_name = tab_file.name.replace('.txt', '.bed')

        # if the output file exists delete it and proceed
        try:
            os.remove(output_name)
        except OSError:
            pass

        working_batch = []
        # assigns the results of get_batch to 'batch' as long as there are lines left to read
        while batch := get_batch(reader, file):
            # keep at least nine lines in tile
            working_batch.extend(batch)
               
            # call the write bed function
            write_bed_file(output_name, working_batch, args)


def get_batch(reader, file):
    return [row for _ in range(9) if (row:=next(reader, None))]


def write_bed_file(output_name, working_batch, args):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chr', 'start', 'end', 'strand', 'coverage1', 'freqC1', 'coverage2', 'freqC2', 'coverage3', 'freqC3']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        # write header to output file if the header flag is true
        if args.header:
            writer.writeheader()

        position = [int(row.get('base')) for row in working_batch]

        while len(working_batch) >= 3:  # this ma cuase the last location to be sckipped if it consists of a range of 2
            lines = []

            if position[0] + 1 != position[1]:
                working_batch.pop(0)
                position.pop(0)
            elif position[0] + 1 == position[1] and position[1] + 1 == position[2]:
                for _ in range(3):
                    lines.append(working_batch.pop(0))
                    position.pop(0)
            elif position[0] + 1 == position[1] and position[1] + 1 != position[2]:
                for _ in range(2):
                    lines.append(working_batch.pop(0))
                    position.pop(0)
            
            if lines:
                line = transform_line(lines)

                # write modified line to the output file in append
                writer.writerow(line)


def transform_line(lines):
    # take the first line to create our transformed line
    line = lines[0]
    length = len(lines)

    # a list of columns from the input file to be removed
    keys_to_remove = ['chrBase', 'freqT']

    # remove undesired columns from read in line
    for k in keys_to_remove:
        line.pop(k, None)

    # rename and merge the dictionary items together into 'line'
    line['start'] = line.pop('base')  # rename column 'base' to 'start'
    line['end'] = int(lines[length-1].get('base')) + 1

    # covert the 'F' and 'R' strand notation for '+' and '-'
    if line.get('strand') == 'F':
        line['strand'] = '+'
    else:
        line['strand'] = '-'

    # add all 'coverage' and 'freqC' column information to the one line
    line['coverage1'] = line.pop('coverage')
    line['freqC1'] = line.pop('freqC')
    line['coverage2'] = lines[1].pop('coverage')
    line['freqC2'] = lines[1].pop('freqC')
    
    # if the length of the genomic range is 3 then retain the methylation info
    if length == 3:
        line['coverage3'] = lines[length-1].pop('coverage')
        line['freqC3'] = lines[length-1].pop('freqC')
    else:
        line['coverage3'] = 0
        line['freqC3'] = 0

    return line


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
        txt_to_bed_file(file, args)


if __name__ == "__main__":
    main()
