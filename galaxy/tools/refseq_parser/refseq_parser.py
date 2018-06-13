#!/usr/bin/python

# Usage: python parser.py [OPTIONS] -r <reference genome> -i <report location>

import sys
import csv
import gzip
import argparse as ap

parser = ap.ArgumentParser(prog='outlier-parser', conflict_handler='resolve',
                           description="Parses output file to eliminate outliers")


input = parser.add_argument_group('Input', '')
input.add_argument('-i', '--input', required=True, metavar="STRING", help="Reports")
input.add_argument('-r', '--reference', required=True, metavar="STRING", help="Reference genome in FASTA format")

if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

"""Determining organism"""
try:
    meta = gzip.open(args.reference, 'rb').readline().split()
except:
    meta = open(args.reference, 'rb').readline().split()
organism = "%s %s" % (meta[1].decode('utf-8'), meta[2].decode('utf-8'))

outlier_flag = False
accession = ""

with open(args.input) as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    next(reader, None)

    for row in reader:
        if row['sample'] == accession and outlier_flag:
            continue
        accession = row['sample']
        taxonomy_split = row['top_taxonomy_name'].split()
        taxonomy = "%s %s" % (taxonomy_split[0], taxonomy_split[1])

        if taxonomy.lower() != organism.lower() or float(row['distance']) > 1e-3:
            print(accession)
            outlier_flag = True