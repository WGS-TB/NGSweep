#!/usr/bin/python

#############################################################
# Dataset preprocessing pipeline                            #
#                                                           #
# Dependencies:                                             #
# Refseq_masher                                             #
# Multiqc                                                   #
# Samtools                                                  #
# Trim_galore                                               #
# Smalt/BWA/Bowtie                                          #
# Kraken                                                    #
#                                                           #
# Ensure all dependencies can be called from the command    #
# line (ex. which)                                          #
#############################################################


import sys
import os
import glob
import subprocess
import logging
import argparse as ap
import preprocess
from mpi4py import MPI


"""Command line interface"""
if __name__ == '__main__':
    # MPI setup
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    parser = ap.ArgumentParser(prog='preprocessing-pipeline', conflict_handler='resolve',
                               description="Preprocessing pipeline - Eliminate outliers from datasets")

    modes = parser.add_argument_group('Modes', '')
    modes.add_argument('--outlier', action='store_true', help="Check for outliers")
    modes.add_argument('--trim', action='store_true', help="Trim FASTQ files")
    modes.add_argument('--map', action='store_true', help="Map FASTQ files to reference and verify with Qualimap")
    modes.add_argument('--kraken', action='store_true', help="Run Kraken and trim contaminated reads")

    kraken = parser.add_argument_group('Kraken Parameters', '')
    kraken.add_argument('--db', metavar="STRING", help="Location of Kraken Database")
    kraken.add_argument('--taxon_id', metavar="STRING", help="Defined taxonomic ID (Default: Match taxon of reference)")

    input = parser.add_argument_group('Input', '')
    input.add_argument('-i', '--input', required=True, metavar="STRING", help="Input dataset (location of all FASTQ files)")
    input.add_argument('-r', '--reference', required=True, metavar="STRING", help="Reference genome in FASTA format")
    input.add_argument('-p', '--paired', action='store_true', help="Paired end reads")

    output = parser.add_argument_group('Output', '')
    output.add_argument('-o', '--outdir', metavar="STRING", help="Output directory containing all reports and outliers")
    output.add_argument('-l', '--log', action='store_true', help="Output a log file")
    output.add_argument('--keepfiles', action='store_true', help="Keep intermediate files")

    aligners = parser.add_argument_group('Aligner', '')
    aligners.add_argument('--bwa', action='store_true', help="Run BWA aligner [Default]")
    aligners.add_argument('--smalt', action='store_true', help="Run Smalt aligner")

    optional = parser.add_argument_group('Optional', '')
    optional.add_argument('-v', '--verbose', action='store_true', help="Print status updates during run to stdout"
                                                                       "use -l to keep a log file instead")
    optional.add_argument('-h','--help', action='help', help="Show help message")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()
    error = 0

    # If no output directory was specified, create an output directory
    if not args.outdir:
        args.outdir = os.path.join(os.getcwd(), "output")

    if not os.path.exists(args.outdir):
        subprocess.call(['mkdir', '-p', args.outdir])

    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        print("The output path is not a directory.")
        sys.exit(2)

    # stdout setup
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")

    stdout = logging.StreamHandler()
    stdout.setFormatter(formatter)
    logger.addHandler(stdout)

    # Log file setup
    if args.log:
        logger.info("Log file initiated")
        log = logging.FileHandler(os.path.join(args.outdir, 'log%d.txt' % rank))
        log.setFormatter(formatter)
        logger.addHandler(log)
        args.verbose = True

    # Check arguments
    if not os.path.isdir(args.input):
        error += 1
        logger.error("Input directory %s does not exist." % args.input)

    if not os.path.isfile(args.reference):
        error += 1
        logger.error("Reference at %s does not exist." % args.reference)

    if not args.outlier and not args.trim and not args.map:
        error += 1
        logger.error("Please specify a mode.")

    if error:
        logger.debug("\nUse --help for more information.")
        parser.print_usage()
        sys.exit(2)

    """Run pipeline"""
    directory = os.path.abspath(args.input)

    # Create directories
    if args.outlier:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'outliers')])
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'mash')])
    if args.trim:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'trimmed_fastq')])
    if args.map:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'mapping')])
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'qualimap')])
    if args.kraken:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'kraken')])
    subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'reports')])

    # Determining goal organism
    organism = preprocess.parseReference(args.verbose, args.reference)

    if rank == 0:
        # Index reference using smalt
        if args.map:
            if args.smalt:
                preprocess.smalt_index(args.reference, args.verbose, args.outdir)
            else:
                preprocess.bwa_index(args.reference, args.verbose, args.outdir)

        jobs = [file for file in os.listdir(directory)]
        jobs = [jobs[i::size] for i in range(size)]
    else:
        jobs = None

    jobs = comm.scatter(jobs, root=0)

    # for file in os.listdir(directory):
    for file in jobs:
        filename = file

        if args.paired and "_2" in filename:
            continue

        if args.paired:
            underscore_position = filename.find('_')

            if underscore_position == -1: # The file is not a paired file
                continue

            accession = filename[:underscore_position]
            dot_position = filename.find('.')
            file_type = filename[dot_position:]

        else:
            dot_position = filename.find('.')
            file_type = filename[dot_position:]
            accession = filename[:dot_position]

        if file_type.lower() not in ['.fasta', '.fastq', '.fasta.gz', '.fastq.gz', '.fq.gz']:
            continue

        if args.verbose:
            logger.info("Working on %s" % accession)
        
        if args.paired:
            if preprocess.check_pairs(accession, file_type, args.input, args.verbose) == 1:
                continue

            preprocess.check_headers(accession, file_type, args.input, args.verbose)

            pipeline = preprocess.preprocess(organism, os.path.join(directory, accession+'_1'+file_type), accession, args.outdir,
                                  args.reference, args.paired, os.path.join(directory, accession+'_2'+file_type),
                                  args.log, args.verbose, args.map, args.outlier, args.trim, args.kraken, args.db,
                                  args.taxon_id)

        else:
            pipeline = preprocess.preprocess(organism, os.path.join(directory, filename), accession, args.outdir,
                                  args.reference, args.paired, None, args.log, args.verbose, args.map, args.outlier,
                                  args.trim, args.kraken, args.db, args.taxon_id)

        if args.outlier:
            pipeline.refseq_masher()
            pipeline.parser(refseq=True, qualimap=False)

            # Check if the file was an outlier
            if glob.glob(os.path.join(os.path.join(args.outdir, 'outliers'), accession)+'*'):
                continue

        if args.trim:
           pipeline.trim_galore()

        if args.kraken:
            pipeline.run_kraken()

        if args.map:
            if args.smalt:
                pipeline.smalt_map()
            else:
                pipeline.bwa_map()
            pipeline.samtools()
            pipeline.qualimap()
            pipeline.parser(refseq=False, qualimap=True)

    preprocess.multiqc(args.verbose, args.outdir)

    if not args.keepfiles:
        preprocess.cleanup(args.verbose, args.outdir)
