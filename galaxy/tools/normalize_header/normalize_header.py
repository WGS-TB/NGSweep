#!/usr/bin/python

# Usage: python parser.py [OPTIONS] -i [reports]

import sys
import re
import os
import gzip
import argparse as ap

"""Check mate pair headers (in case fastq files were not downloaded with the -origfmt flag using fastq-dump)"""
def check_headers(accession, file_type, input):
    firstLineForward = gzip.open(os.path.join(input, accession+"_1"+file_type), 'r').readline()
    firstLineBackward = gzip.open(os.path.join(input, accession+"_2"+file_type), 'r').readline()

    if firstLineForward == firstLineBackward: # If the headers are the same (unmodified)
        return

    forwardPaired = firstLineForward.find("/1".encode())
    backwardPaired = firstLineBackward.find("/2".encode())

    if firstLineForward[:forwardPaired] == firstLineBackward[:backwardPaired]: # Paired read identifiers present
        return

    if re.search('(\.1){2}', firstLineForward.decode('utf-8')):
        for i in range(1,3):
            new_file = ""
            with gzip.open(os.path.join(input, accession+"_%d"+file_type) % i, 'r') as mate:
                for line in mate.readlines():
                    line = line.decode('utf-8')
                    if "%s" % accession in line:
                        split = line.split()
                        header = "%s %s" % (split[0][:-2], split[1]+"/%d\n" % i)
                        new_file += header
                        continue
                    new_file += line
            with gzip.open(os.path.join(input, accession+"_%d"+file_type) % i, 'wb') as mate:
                mate.write(str.encode(new_file))


parser = ap.ArgumentParser(prog='outlier-parser', conflict_handler='resolve',
                           description="Parses output file to eliminate outliers")

input = parser.add_argument_group('Input', '')
input.add_argument('-i', '--input', nargs='+', required=True, help="Tab-delimited RefSeq Masher reports")

if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

if not os.path.exists('output'):
    os.makedirs('output')

for file in args.input:
    filename = os.path.split(file)[1]

    underscore_position = filename.find('_')
    accession = filename[:underscore_position]
    dot_position = filename.find('.')
    file_type = filename[dot_position:]

    check_headers(accession, file_type, file)

