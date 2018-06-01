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
import gzip
import csv
import logging
import argparse as ap


class preprocess():

    def __init__(self, organism, input, name, outdir, reference, paired, input2, log, verbose, map, outlier, trim,
                 kraken, db, taxon_id):
        self.organism = organism
        self.input = input
        self.name = name
        self.outdir = outdir
        self.reference = reference
        self.paired = paired
        self.input2 = input2
        self.log = log
        self.verbose = verbose
        self.outlier = outlier
        self.map = map
        self.trim = trim
        self.kraken = kraken
        self.db = db
        self.taxon_id = taxon_id
        self.logger = logging.getLogger()

    """Shell Execution"""
    def runCommand(self, command, directory):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=directory)
        out,err = process.communicate()

        if out:
            self.logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
        if err:
            self.logger.info("Standard error: \n" + err.decode('utf-8') + "\n")

    """Running Refseq_masher"""
    def refseq_masher(self):
        self.ifVerbose("Running Refseq_masher matching")
        self.runCommand(['refseq_masher', 'matches', '-o', self.name+'.match', '--output-type', 'tab', self.input],
                        os.path.join(self.outdir, 'mash'))

    """Running Kraken"""
    def run_kraken(self):
        self.ifVerbose("Running Kraken")

        if self.paired:
            self.runCommand(['kraken', '--db', self.db, '--classified-out', self.name + '.classified', '--paired',
                             '--output', self.name + '.kraken', self.input, self.input2],
                             os.path.join(args.outdir, 'kraken'))
        else:
            self.runCommand(['kraken', '--db', self.db, '--classified-out', self.name + '.classified', '--output',
                             self.name + '.kraken', self.input], os.path.join(args.outdir, 'kraken'))

    """Run Trim_galore to preprocess fastq files"""
    def trim_galore(self):
        self.ifVerbose("Trimming fastq files using Trim_galore")

        if self.paired:
            self.runCommand(['trim_galore', '--fastqc_args', "\"--outdir " +
                             os.path.join(self.outdir, "trimmed_fastq/fastqc") + "\"", '--gzip', '-o',
                             os.path.join(self.outdir, "trimmed_fastq"), '--paired', self.input, self.input2], None)

            self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_1_val_1.fq.gz")
            self.input2 = self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_2_val_2.fq.gz")

        else:
            self.runCommand(['trim_galore', '--fastqc', '--gzip', '-o',
                             os.path.join(self.outdir, "trimmed_fastq"), self.input], None)

            self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_val.fq.gz")

    """Mapping with Smalt"""
    def smalt_map(self):
        self.ifVerbose("Mapping reads to reference")
        if self.paired:
            self.runCommand(['smalt', 'map', '-i', '1000', '-j', '20', '-l', 'pe', '-o', self.name+".BAM",
                             'reference', self.input, self.input2], os.path.join(self.outdir, 'bam'))
        else:
            self.runCommand(['smalt', 'map', '-o', self.name+".BAM", 'reference', self.input],
                            os.path.join(args.outdir, 'bam'))

    """Sort BAM files using Samtools"""
    def samtools(self):
        self.ifVerbose("Sorting BAM files using Samtools")
        self.runCommand(['samtools', 'sort', '-o', self.name+'_sorted.BAM', self.name+'.BAM'],
                        os.path.join(args.outdir, 'bam'))

    """Checking mapping quality with Qualimap"""
    def qualimap(self):
        self.ifVerbose("Running qualimap BAM QC")
        self.runCommand(['qualimap', 'bamqc', '-bam', 'temp/'+self.name+'_sorted.BAM', '-outformat', 'HTML',
                         '-outdir ', os.path.join(self.outdir, self.name)], None)

    """Parse through report file obtained from Qualimap or Refseq_masher"""
    def parser(self):
        if self.outlier:
            outlier_flag = False

            self.ifVerbose("Parsing Refseq_masher report")
            with open(os.path.join(os.path.join(args.outdir, 'mash'), self.name + '.match')) as csvfile:
                for row in csv.DictReader(csvfile, delimiter='\t'):
                    taxonomy_split = row['top_taxonomy_name'].split()
                    taxonomy = "%s %s" % (taxonomy_split[0], taxonomy_split[1])

                    if taxonomy.lower() != self.organism.lower() or float(row['distance']) > 1e-3:
                        outlier_flag = True
                        break

            if self.map and not outlier_flag:
                self.ifVerbose("Parsing Qualimap report")
                report = open(self.outdir+'/'+self.name+'/genome_results.txt')

                for line in report:
                    if "number of mapped reads" in line:
                        mapped_percentage = line.split()[-1].strip('()%')
                    if "mean mapping quality" in line:
                        mean_mapping_quality = line.split()[-1]

                if float(mapped_percentage) < 90 or float(mean_mapping_quality) < 10:
                    outlier_flag = True

            if outlier_flag:
                self.ifVerbose("%s is an outlier" % self.name)
                self.move_outlier()

    """Move outlier sample"""
    def move_outlier(self):
        self.ifVerbose("Moving %s to outlier folder" % self.name)
        self.runCommand(['cp', self.input, os.path.join(self.outdir, 'outliers')], None)
        if self.paired:
            self.runCommand(['cp', self.input2, os.path.join(self.outdir, 'outliers')], None)

    def ifVerbose(self, msg):
        if self.verbose:
            self.logger.info(msg)

"""Check mate pair consistency"""
def check_pairs(accession, file_type, input, verbose):
    logger = logging.getLogger()
    if verbose:
        logger.info("Checking consistency of mate pairs")

    files = glob.glob(os.path.join(input, accession + "*"))
    if os.path.join(input, accession+"_1"+file_type) not in files or os.path.join(input, accession+"_2"+file_type) not in files:
        logger.error("%s is missing a mated pair" % accession)
        return 1
    logger.info("OK")
    return 0

"""Indexing reference with Smalt"""
def smalt_index(reference, verbose, outdir):
    logger = logging.getLogger()

    if verbose:
        logger.info("Creating index file for reference sequence")

    process = subprocess.Popen(['smalt', 'index', '-k', '13', '-s', '6', 'reference', reference],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.join(outdir, 'bam'))
    out,err = process.communicate()

    if out:
        logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
    if err:
        logger.info("Standard error: \n" + err.decode('utf-8') + "\n")

"""Determine organism from reference"""
def parseReference(verbose, reference):
    logger = logging.getLogger()

    if verbose:
        logger.info("Determining organism from reference given")
    meta = gzip.open(reference, 'rb').readline().split()
    organism = "%s %s" % (meta[1].decode('utf-8'), meta[2].decode('utf-8'))
    return organism

"""Generating reports with MultiQC"""
def multiqc(verbose, outdir):
    logger = logging.getLogger()

    if verbose:
        logger.info("Generating reports using MultiQC")

    process = subprocess.Popen(['multiqc', '.', '-o', os.path.join(outdir, 'reports')],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=outdir)
    out, err = process.communicate()

    if out:
        logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
    if err:
        logger.info("Standard error: \n" + err.decode('utf-8') + "\n")

"""Removes all temporary files"""
def cleanup(verbose, outdir):
    logger = logging.getLogger()

    if verbose:
        logger.info("Deleting temporary files")

    temp_folders = ['mash', 'bam', 'qualimap']

    for folder in temp_folders:
        process = subprocess.Popen(['rm', '-r', folder],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=outdir)
        out, err = process.communicate()

        if out:
            logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
        if err:
            logger.info("Standard error: \n" + err.decode('utf-8') + "\n")


"""Command line interface"""
if __name__ == '__main__':
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
    aligners.add_argument('--smalt', action='store_true', help="Run the Smalt aligner (Default)")
    # aligners.add_argument('--bowtie', action='store_true', help="Run Bowtie aligner")

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
        log = logging.FileHandler(os.path.join(args.outdir, 'log.txt'))
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

    # Todo: Give a choice of aligners

    """Run pipeline"""
    directory = os.path.abspath(args.input)

    # Create directories
    if args.outlier:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'outliers')])
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'mash')])
    if args.trim:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'trimmed_fastq')])
    if args.map:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'bam')])
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'qualimap')])
    if args.kraken:
        subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'kraken')])
    subprocess.call(['mkdir', '-p', os.path.join(args.outdir, 'reports')])

    # Index reference using smalt
    if args.map:
        smalt_index(args.reference, args.verbose, args.outdir)

    # Determining goal organism
    organism = parseReference(args.verbose, args.reference)

    for file in os.listdir(directory):
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
            if check_pairs(accession, file_type, args.input, args.verbose) == 1:
                continue

            pipeline = preprocess(organism, os.path.join(directory, accession+'_1'+file_type), accession, args.outdir,
                                  args.reference, args.paired, os.path.join(directory, accession+'_2'+file_type),
                                  args.log, args.verbose, args.map, args.outlier, args.trim, args.kraken, args.db,
                                  args.taxon_id)

        else:
            pipeline = preprocess(organism, os.path.join(directory, filename), accession, args.outdir,
                                  args.reference, args.paired, None, args.log, args.verbose, args.map, args.outlier,
                                  args.trim, args.kraken, args.db, args.taxon_id)

        if args.outlier:
            pipeline.refseq_masher()
            pipeline.parser()

            # Check if the file was an outlier
            if glob.glob(os.path.join(os.path.join(args.outdir, 'outliers'), accession)+'*'):
                continue

        if args.trim:
           pipeline.trim_galore()

        if args.kraken:
            pipeline.run_kraken()

        if args.map:
            pipeline.smalt_map()
            pipeline.samtools()
            pipeline.qualimap()

    multiqc(args.verbose, args.outdir)

    if not args.keepfiles:
        cleanup(args.verbose, args.outdir)
