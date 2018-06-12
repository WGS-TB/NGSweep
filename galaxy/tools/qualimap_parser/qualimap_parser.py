#!/usr/bin/python

# Usage: python parser.py [OPTIONS] -i <report location>

import sys
import argparse as ap

parser = ap.ArgumentParser(prog='outlier-parser', description="Parses output file to eliminate outliers")


input = parser.add_argument_group('Input', '')
input.add_argument('-i', '--input', required=True, metavar="STRING", help="Reports")

if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

"""Open output file"""
output = open("outlier_list_qualimap.csv", "w")

report = open(args.input)

for line in report:
    if "number of mapped reads" in line:
        mapped_percentage = line.split()[-1].strip('()%')
    if "mean mapping quality" in line:
        mean_mapping_quality = line.split()[-1]

    if float(mapped_percentage) < 90 or float(mean_mapping_quality) < 10:
        outlier_flag = True