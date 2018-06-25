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