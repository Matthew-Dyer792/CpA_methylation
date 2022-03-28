# file stream script to analyze CpA methylation
from asyncore import read
import csv
import os

with open('DNMT_practice.txt', 'r', buffering=5368709120) as tab_file:  # 1GB = 1073741824 5GB = 5368709120
    reader = csv.DictReader(tab_file, delimiter='\t')

    output_name = f"{tab_file.name.split('.')[0]}.bed"

    try:
        os.remove(output_name)
    except OSError:
        pass


    with open(output_name, 'a') as output_file:
        field_names = ['chr', 'start', 'end', 'strand']

        writer = csv.DictWriter(output_file, fieldnames=field_names, delimiter='\t')

        writer.writeheader()

        for line in reader:
            keys_to_remove = ['chrBase', 'coverage', 'freqC', 'freqT']

            for k in keys_to_remove:
                line.pop(k, None)

            line['start'] = line.pop('base')

            line['end'] = int(line.get('start')) + 1

            if line.get('strand') == 'F':
                line['strand'] = '+'
            else:
                line['strand'] = '-'

            writer.writerow(line)
            # print(line)
