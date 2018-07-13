#!/usr/bin/python

import os
import subprocess
import gzip
import glob
import csv
import logging
import re
from ete3 import NCBITaxa

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
    def runCommand(self, command, directory, write_output):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=directory)
        out,err = process.communicate()

        if out:
            if write_output:
                return out
            self.logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
        if err:
            self.logger.info("Standard error: \n" + err.decode('utf-8') + "\n")

    """Running Refseq_masher"""
    def refseq_masher(self):
        self.ifVerbose("Running Refseq_masher matching")
        self.runCommand(['refseq_masher', 'matches', '-o', self.name+'.match', '--output-type', 'tab', self.input],
                        os.path.join(self.outdir, 'mash'), write_output=False)

    """Running Kraken"""
    def run_kraken(self):
        self.ifVerbose("Running Kraken")

        gzip = ""
        if self.input[-3:] == ".gz":
            gzip = "--gzip-compressed"

        if self.paired:
            self.runCommand(['kraken', '--db', self.db, '--paired', '--output', self.name + '.kraken',
                             '--fastq-input', "%s" % gzip, self.input, self.input2],
                             os.path.join(self.outdir, 'kraken'), write_output=False)
        else:
            self.runCommand(['kraken', '--db', self.db, '--output', self.name + '.kraken', '--fastq-input',
                             "%s" % gzip, self.input], os.path.join(self.outdir, 'kraken'), write_output=False)

    """Parse Kraken resuts"""
    def parse_kraken_results(self):
        self.ifVerbose("Parsing Kraken results")

        kraken = {} # Store classification for each read

        with open(os.path.join(self.outdir, 'kraken/%s.kraken' % self.name), 'r') as classification:
            for line in classification:
                classified, read_id, tax_id, length, details = line.strip().split("\t")
                kraken[read_id] = tax_id

        # Obtain taxonomic ID for descendants of target organism
        ncbi = NCBITaxa()
        descendants = ncbi.get_descendant_taxa(self.taxon_id)

        # Classify each read
        kraken_class = {}

        for read_id, tax_id in kraken.iteritems():
            if tax_id == 0:
                kraken_class[read_id] = "unclassified"
            elif tax_id in descendants:
                kraken_class[read_id] = "target"
            else:
                kraken_class[read_id] = "other"

        return kraken_class

    """Trim fastq reads not belonging to target organism"""
    def kraken_trim(self):
        self.ifVerbose("Trimming reads that do not belong to the target organism")

        kraken = self.parse_kraken_results()

        # Write new fastq file
        if self.paired:
            files = glob.glob(os.path.join(self.input, self.name + "*"))
            for fastq_in in files:
                with gzip.open(fastq_in) as f_in:
                    fastq_out = os.path.split(file)[1]
                    if fastq_out[-3:] == ".gz": # Eliminate .gz from filename
                        fastq_out = fastq_out[:-3]
                    with open(os.path.join(self.outdir, 'kraken_trim/%s' % fastq_out)) as f_out:
                        for line in f_in:
                            # Split ID with space, then remove "/1" or "/2" if it exists and ignore initial @
                            read_id = line.split(" ")[0].split("/")[0][1:]
                            if read_id in kraken and kraken[read_id] != "other":
                                f_out.write(line)
                                for i in range(3):
                                    f_out.write(f_in.readline())
                            else:
                                for i in range(3):
                                    f_in.readline()

                # Zip output files
                self.runCommand(['gzip', os.path.join(self.outdir, 'kraken_trim/%s' % fastq_out)], None, write_output=False)


    """Run Trim_galore to preprocess fastq files"""
    def trim_galore(self):
        self.ifVerbose("Trimming fastq files using Trim_galore")

        if self.paired:
            self.runCommand(['trim_galore', '--fastqc_args', "\"--outdir " +
                             os.path.join(self.outdir, "trimmed_fastq/fastqc") + "\"", '--gzip', '-o',
                             os.path.join(self.outdir, "trimmed_fastq"), '--paired', self.input, self.input2],
                            directory=None, write_output=False)

            self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_1_val_1.fq.gz")
            self.input2 = self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_2_val_2.fq.gz")

        else:
            self.runCommand(['trim_galore', '--fastqc', '--gzip', '-o',
                             os.path.join(self.outdir, "trimmed_fastq"), self.input], directory=None, write_output=False)

            self.input = os.path.join(os.path.join(self.outdir, "trimmed_fastq"), self.name + "_val.fq.gz")

    """Mapping with Smalt"""
    def smalt_map(self):
        self.ifVerbose("Mapping reads to reference using Smalt")
        if self.paired:
            self.runCommand(['smalt', 'map', '-i', '1000', '-j', '20', '-l', 'pe', '-o', self.name+".BAM",
                             'reference', self.input, self.input2], os.path.join(self.outdir, 'mapping'), write_output=False)
        else:
            self.runCommand(['smalt', 'map', '-o', self.name+".BAM", 'reference', self.input],
                            os.path.join(self.outdir, 'mapping'), write_output=False)

    """Mapping with BWA"""
    def bwa_map(self):
        self.ifVerbose("Mapping reads to reference using BWA")
        with open(os.path.join(self.outdir, 'mapping/%s.SAM' % self.name), 'wb') as sam:
            with open(os.path.join(self.outdir, 'mapping/%s.BAM' % self.name), 'wb') as bam:
                if self.paired:
                    sam_output = self.runCommand(['bwa', 'mem', 'reference', self.input, self.input2],
                                                 os.path.join(self.outdir, 'mapping'), write_output=True)
                    sam.write(sam_output)
                    bam_output = self.runCommand(['samtools', 'view', '-Sb', self.name+".SAM"],
                                                 os.path.join(self.outdir, 'mapping'), write_output=True)
                    bam.write(bam_output)
                else:
                    sam_output = self.runCommand(['bwa', 'mem', 'reference', self.input, '>', self.name+".SAM"],
                                                 os.path.join(self.outdir, 'mapping'), write_output=True)
                    sam.write(sam_output)
                    bam_output = self.runCommand(['samtools', 'view', '-Sb', self.name+".SAM", '>', self.name+".BAM"],
                                                 os.path.join(self.outdir, 'mapping'), write_output=True)
                    bam.write(bam_output)

    """Sort BAM files using Samtools"""
    def samtools(self):
        self.ifVerbose("Sorting BAM files using Samtools")
        self.runCommand(['samtools', 'sort', '-o', self.name+'_sorted.BAM', self.name+'.BAM'],
                        os.path.join(self.outdir, 'mapping'), write_output=False)

    """Checking mapping quality with Qualimap"""
    def qualimap(self):
        self.ifVerbose("Running qualimap BAM QC")

        self.runCommand(['qualimap', 'bamqc', '-bam', os.path.join(self.outdir, 'mapping/'+self.name+'_sorted.BAM'), '-outformat', 'HTML',
                         '-outdir ', os.path.join(self.outdir, "qualimap/"+self.name)], directory=None, write_output=False)


    """Parse through report file obtained from Qualimap or Refseq_masher"""
    def parser(self, refseq, qualimap):
        if self.outlier:
            outlier_flag = False

            if refseq:
                self.ifVerbose("Parsing Refseq_masher report")
                with open(os.path.join(os.path.join(self.outdir, 'mash'), self.name + '.match')) as csvfile:
                    for row in csv.DictReader(csvfile, delimiter='\t'):
                        taxonomy_split = row['top_taxonomy_name'].split()
                        taxonomy = "%s %s" % (taxonomy_split[0], taxonomy_split[1])

                        if taxonomy.lower() != self.organism.lower() or float(row['distance']) > 1e-3:
                            outlier_flag = True
                            break

            if qualimap:
                self.ifVerbose("Parsing Qualimap report")
                report = open(self.outdir+'/qualimap/'+self.name+'/genome_results.txt')

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
        self.runCommand(['cp', self.input, os.path.join(self.outdir, 'outliers')], directory=None, write_output=False)
        if self.paired:
            self.runCommand(['cp', self.input2, os.path.join(self.outdir, 'outliers')], directory=None, write_output=False)

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

"""Check mate pair headers (in case fastq files were not downloaded with the -origfmt flag using fastq-dump)"""
def check_headers(accession, file_type, input, verbose):
    logger = logging.getLogger()

    if verbose:
        logger.info("Checking mate pair headers")

    flag = False

    with gzip.open(os.path.join(input, accession+"_1"+file_type), 'r') as forward:
        firstLineForward = forward.readline()
    with gzip.open(os.path.join(input, accession+"_2"+file_type), 'r') as reverse:
        firstLineBackward = reverse.readline()

    if firstLineForward == firstLineBackward: # If the headers are the same (unmodified)
        return 0

    forwardPaired = firstLineForward.find("/1".encode())
    backwardPaired = firstLineBackward.find("/2".encode())

    if firstLineForward[:forwardPaired] == firstLineBackward[:backwardPaired]: # Paired read identifiers present
        return 0

    if re.search('(\.1){2}', firstLineForward.decode('utf-8')):
        if flag and verbose:
            logger.info("The headers are modified NCBI SRA identifiers")

        for i in range(1,3):
            new_file = ""
            with gzip.open(os.path.join(input, accession+"_%d"+file_type) % i, 'r') as mate:
                for line in mate.readlines():
                    line = line.decode('utf-8')
                    if "%s" % accession in line:
                        flag = True
                        split = line.split()
                        header = "%s %s" % (split[0][:-2], split[1]+"/%d\n" % i)
                        new_file += header
                        continue
                    new_file += line
            with gzip.open(os.path.join(input, accession+"_%d"+file_type) % i, 'wb') as mate:
                mate.write(str.encode(new_file))

        if flag and verbose:
            logger.info("Headers have been edited to ensure uniformity")
        return 0

    return 1


"""Indexing reference with Smalt"""
def smalt_index(reference, verbose, outdir):
    logger = logging.getLogger()

    if verbose:
        logger.info("Creating index file for reference sequence using Smalt")

    process = subprocess.Popen(['smalt', 'index', '-k', '13', '-s', '6', 'reference', reference],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.join(outdir, 'mapping'))
    out,err = process.communicate()

    if out:
        logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
    if err:
        logger.info("Standard error: \n" + err.decode('utf-8') + "\n")

def bwa_index(reference, verbose, outdir):
    logger = logging.getLogger()

    if verbose:
        logger.info("Creating index file for reference sequence using BWA")

    process = subprocess.Popen(['bwa', 'index', '-p', 'reference', reference],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.join(outdir, 'mapping'))
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

    temp_folders = ['mash', 'mapping', 'qualimap']

    for folder in temp_folders:
        process = subprocess.Popen(['rm', '-r', folder],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=outdir)
        out, err = process.communicate()

        if out:
            logger.info("Standard output: \n" + out.decode('utf-8') + "\n")
        if err:
            logger.info("Standard error: \n" + err.decode('utf-8') + "\n")
