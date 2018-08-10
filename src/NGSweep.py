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
# ete3                                                      #
#                                                           #
# Ensure all dependencies can be called from the command    #
# line (ex. which)                                          #
#############################################################


import sys
import os
import glob
import logging
import argparse as ap
import subprocess as sp
import preprocess
from mpi4py import MPI

VERSION = "0.1"

"""Create directory"""
def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

"""Command line interface"""
if __name__ == '__main__':
    parser = ap.ArgumentParser(prog='ngsweep', conflict_handler='resolve',
                               usage='%(prog)s [OPTIONS] {outlier,trim,kraken,qa}',
                               description="NGSweep: Next-Generation Sequencing Data Preprocessing Pipeline\n"
                                           "Version %s\n"
                                           "Contact: Matthew Nguyen <mtn14@sfu.ca>" % VERSION,
                               formatter_class=ap.RawTextHelpFormatter)

    parent_parser = ap.ArgumentParser(add_help=False)
    parent_parser.add_argument('inputdir', help='Input dataset (location of all FASTQ files)')
    parent_parser.add_argument('-p', '--paired', action='store_true', help="Paired end reads")

    parent_parser.add_argument('-o', '--outdir', metavar="DIR",
                               help="Output directory containing all reports and outliers")
    parent_parser.add_argument('-l', '--log', action='store_true', help="Output a log file")
    parent_parser.add_argument('-kf', '--keepfiles', action='store_true', help="Keep intermediate files")
    parent_parser.add_argument('-v', '--verbose', action='store_true', help="Print status updates to stdout")

    subparser = parser.add_subparsers(help="commands")

    outlier_parser = subparser.add_parser('outlier', parents=[parent_parser],
                                          usage='ngsweep outlier [OPTIONS] <inputdir> <taxid>',
                                          help='Find outliers in an NGS dataset')
    outlier_parser.add_argument('taxid', type=int, help="Target taxonomic ID")

    trim_parser = subparser.add_parser('trim', parents=[parent_parser],
                                       usage='ngsweep trim [OPTIONS]',
                                       help="Trim adapters and low quality bases from FASTQ sequences")

    kraken_parser = subparser.add_parser('kraken', parents=[parent_parser],
                                         usage='ngsweep kraken [OPTIONS] <inputdir> <db> <taxid>',
                                         help="Trim contaminated reads from FASTQ sequences")
    kraken_parser.add_argument('db', help="Location of Kraken Database")
    kraken_parser.add_argument('taxid', type=int,
                               help="Target taxonomic ID")

    qa_parser = subparser.add_parser('qa', parents=[parent_parser],
                                     usage='ngsweep qa [OPTIONS] <inputdir>', help="Perform QA on FASTQ samples")

    # input.add_argument('-r', '--reference', metavar="STRING", help="Reference genome in FASTA format")
    #
    # aligners = parser.add_argument_group('Aligner', '')
    # aligners.add_argument('--bwa', action='store_true', help="Run BWA aligner [Default]")
    # aligners.add_argument('--smalt', action='store_true', help="Run Smalt aligner")

    parser.add_argument('--test', action='store_true', help="Test for dependencies")
    parser.add_argument('-h','--help', action='help', help="Show help message")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()

    # stdout setup
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")

    stdout = logging.StreamHandler()
    stdout.setFormatter(formatter)
    logger.addHandler(stdout)

    if args.test:
        logger.info("Testing for dependencies")
        dependencies = ['refseq_masher', 'multiqc', 'samtools', 'trim_galore', 'kraken']
        for dependency in dependencies:
            sp.call([dependency, '--version'])
        sp.call(['smalt','version'])
        sp.call(['bwa'])
        sys.exit()

    error = 0

    # MPI setup
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # If no output directory was specified, create an output directory
    if not args.outdir:
        args.outdir = os.path.join(os.getcwd(), "output")

    if rank == 0:
        mkdir(args.outdir)

    if os.path.exists(args.outdir) and not os.path.isdir(args.outdir):
        logger.error("The output path is not a directory.")
        sys.exit(2)

    comm.Barrier()

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

    if not args.outlier and not args.trim and not args.map and not args.kraken:
        error += 1
        logger.error("Please specify a mode.")

    if error:
        logger.debug("\nUse --help for more information.")
        parser.print_usage()
        sys.exit(2)

    """Run pipeline"""
    directory = os.path.abspath(args.input)

    if rank == 0:
        # Create directories
        if args.outlier:
            # mkdir(os.path.join(args.outdir, 'outliers'))
            mkdir(os.path.join(args.outdir, 'mash'))
        if args.trim:
            mkdir(os.path.join(args.outdir, 'trimmed_fastq'))
        if args.map:
            mkdir(os.path.join(args.outdir, 'mapping'))
            mkdir(os.path.join(args.outdir, 'qualimap'))
        if args.kraken:
            mkdir(os.path.join(args.outdir, 'kraken'))
            mkdir(os.path.join(args.outdir, 'kraken_trim'))
        mkdir(os.path.join(args.outdir, 'reports'))

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

            if preprocess.check_headers(accession, file_type, args.input, args.verbose) == 1:
                logger.critical("%s has an undefined error in the header.  Please check the headers of the paired reads "
                                "to ensure uniformity." % accession)

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
            pipeline.kraken_trim()

        if args.map:
            if args.smalt:
                pipeline.smalt_map()
            else:
                pipeline.bwa_map()
            pipeline.samtools()
            pipeline.qualimap()
            pipeline.parser(refseq=False, qualimap=True)
            
    comm.Barrier()

    preprocess.multiqc(args.verbose, args.outdir)

    if not args.keepfiles:
        preprocess.cleanup(args.verbose, args.outdir)

    # Conglomerate all log files
    log = open("output/log.txt", 'a')

    for index,file in enumerate(glob.glob(os.path.join(args.outdir, "log*.txt"))):
        log.write("PROCESS #%d\n" % index)
        log.write("-" * 25)
        log.write(file+"\n\n")
        os.remove(file)
