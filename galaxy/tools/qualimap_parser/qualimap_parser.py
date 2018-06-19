#!/usr/bin/python

# Usage: python parser.py [OPTIONS] -i [reports]

import sys
import argparse as ap


### Define global variables
accession = ""
mapped_percentage = ""
mean_mapping_quality = ""

parser = ap.ArgumentParser(prog='outlier-parser', conflict_handler='resolve',
                           description="Parses output file to eliminate outliers")

input = parser.add_argument_group('Input', '')
input.add_argument('-i', '--input', nargs='+', required=True, help="Tab-delimited RefSeq Masher reports")

if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()
output = open("outlier_list.tsv", "a")

for file in args.input:
    with open(file, 'r') as report:
        iterable = iter(report)

        for line in iterable:
            if "BAM file:" in line and not accession:
                line = next(iterable)
                accession = line.split('</td>')[0].strip("<td class=column2>")
                continue

            if "Mapped reads" in line and not mapped_percentage:
                line = next(iterable)
                mapped_percentage = line.split('</td')[0].strip("<td class=column2>").split()[2].strip("%")
                continue

            if "Mean Mapping Quality" in line and not mean_mapping_quality:
                line = next(iterable)
                mean_mapping_quality = line.split('</td>')[0].strip("<td class=column2>")
                break

        if float(mapped_percentage) < 90 or float(mean_mapping_quality) < 10:
            output.write("%s\n" % accession[:accession.find(".")])