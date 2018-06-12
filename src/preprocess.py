#!/usr/bin/python

import os
import subprocess
import gzip
import glob
import csv
import logging

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
                             os.path.join(self.outdir, 'kraken'))
        else:
            self.runCommand(['kraken', '--db', self.db, '--classified-out', self.name + '.classified', '--output',
                             self.name + '.kraken', self.input], os.path.join(self.outdir, 'kraken'))

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
                            os.path.join(self.outdir, 'bam'))

    """Sort BAM files using Samtools"""
    def samtools(self):
        self.ifVerbose("Sorting BAM files using Samtools")
        self.runCommand(['samtools', 'sort', '-o', self.name+'_sorted.BAM', self.name+'.BAM'],
                        os.path.join(self.outdir, 'bam'))

    """Checking mapping quality with Qualimap"""
    def qualimap(self):
        self.ifVerbose("Running qualimap BAM QC")
        self.runCommand(['qualimap', 'bamqc', '-bam', os.path.join(self.outdir, 'bam/'+self.name+'_sorted.BAM'), '-outformat', 'HTML',
                         '-outdir ', os.path.join(self.outdir, self.name)], None)

    """Parse through report file obtained from Qualimap or Refseq_masher"""
    def parser(self):
        if self.outlier:
            outlier_flag = False

            self.ifVerbose("Parsing Refseq_masher report")
            with open(os.path.join(os.path.join(self.outdir, 'mash'), self.name + '.match')) as csvfile:
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