# file stream script to analyze CpA methylation
import csv
import os
import glob
import argparse


# file_line = 0

# chr_file_names = {
#     "chr1": f"{}"
# }


def split_file_on_chr(file):
    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        for line in reader:
            # get the chromosome of the current line
            chromosome = line.get("chr")

            # # if the chromosome exists in the keys then write the line
            # if chromosome in chr_file_names.keys():
            #     # open the proper chromosome file to write the line to
            #     with open(f"temp/{chr_file_names.get(chromosome)}.txt", 'a') as output_file:
            #         writer = csv.DictWriter(output_file, delimiter='\t')

            #         writer.writerow(line)


def txt_to_bed_file(file, args):
    # open txt files in excess of 24GB using a buffer size of 2GB
    with open(file, 'r', buffering=int(2.1e9)) as tab_file:
        reader = csv.DictReader(tab_file, delimiter='\t')

        # design output file name based off the input file (substitute .txt for .bed)
        output_name = tab_file.name.replace('.txt', '.bed')

        # if the output file exists delete it and proceed
        try:
            os.remove(output_name)
        except OSError:
            pass

        working_batch = []
        # assigns the results of get_batch to 'batch' as long as there are lines left to read
        while batch := get_batch(reader):
            # keep at least nine lines in tile
            working_batch.extend(batch)

            # call the 'write_bed' function
            write_bed_file(output_name, working_batch, args)


def get_batch(reader):
    # list comprehension to produce a list of the next nine rows as dictionaries
    return [row for _ in range(9) if (row := next(reader, None))]


def write_bed_file(output_name, working_batch, args):
    # open output file in append mode
    with open(output_name, 'a') as output_file:
        # set the names of the fields to be written to the new file
        field_names = ['chr', 'start', 'end', 'strand', 'coverage1', 'freqC1',
                       'coverage2', 'freqC2', 'coverage3', 'freqC3']

        # create line writer
        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        # write header to output file if the header flag is true
        if args.header:
            writer.writeheader()

        # list comprehension to produce an output of all the genomic locations in our tile 
        position = [int(row.get('base')) for row in working_batch]

        # this might cause the last location to be skipped if it consists of a range of 2
        while len(working_batch) >= 3:
            lines = []

            # if the run of consecutive genomic regions is one
            if position[0] + 1 != position[1]:
                # remove the genomic region of one span
                lines.append(working_batch.pop(0))
                # working_batch.pop(0)
                position.pop(0)
            # if the run of consecutive genomic regions is three
            elif position[0] + 1 == position[1] and position[1] + 1 == position[2]:
                # transfer the three relevant lines on
                for _ in range(3):
                    lines.append(working_batch.pop(0))
                    position.pop(0)
            # if the run of consecutive genomic regions is two
            elif position[0] + 1 == position[1] and position[1] + 1 != position[2]:
                # transfer the two relevant lines on
                for _ in range(2):
                    lines.append(working_batch.pop(0))
                    position.pop(0)

            # if lines is not empty hand it off to be flattened to one line
            if lines:
                line = transform_line(lines)

                # write modified line to the output file in append
                writer.writerow(line)


def transform_line(lines):
    # take the first line to create our transformed line
    line = lines[0]
    length = len(lines)
    end = int(lines[length - 1].get('base')) + 1

    # a list of columns from the input file to be removed
    keys_to_remove = ['chrBase', 'freqT']

    # remove undesired columns from read in line
    for k in keys_to_remove:
        line.pop(k, None)

    # rename and merge the dictionary items together into 'line'
    line['start'] = line.pop('base')  # rename column 'base' to 'start'
    line['end'] = end

    # covert the 'F' and 'R' strand notation for '+' and '-'
    if line.get('strand') == 'F':
        line['strand'] = '+'
    else:
        line['strand'] = '-'

    # add all 'coverage' and 'freqC' column information to the one line
    line['coverage1'] = line.pop('coverage')
    line['freqC1'] = line.pop('freqC')

    # if the length of the genomic range is 2 then retain the methylation info
    if length >= 2:
        line['coverage2'] = lines[1].pop('coverage')
        line['freqC2'] = lines[1].pop('freqC')
    else:  # prevent unequal row values
        line['coverage2'] = 0
        line['freqC2'] = 0

    # if the length of the genomic range is 3 then retain the methylation info
    if length == 3:
        line['coverage3'] = lines[length - 1].pop('coverage')
        line['freqC3'] = lines[length - 1].pop('freqC')
    else:  # prevent unequal row values
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
        # split_file_on_chr(file)
        # txt_to_bed_file(file, args)
        print(file)


if __name__ == "__main__":
    main()
