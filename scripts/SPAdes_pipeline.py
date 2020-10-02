#!/usr/bin/env python3

import os
import subprocess
import argparse
import time

"""Notes for an eventual protocol for this pipeline"""
#A module file serving as a config to set up pathways will be included on the server

"""Variables for database pathways; these need to be set until pipeline upgrade"""
adapter_database1 = '/data/opt/jarfiles/trimmomatic/v0.39/Trimmomatic-0.39/adapters/All_adapters.fa'
kraken2_database1 = '/data/opt/databases/kraken2/RefSeq_Kraken2_std_2019_02_12'


"""Full pathway for trimmomatic jar and removesmalls.pl scripts. Not sure how to include in sysargs argument at this point"""

trimmomatic_path = '/data/opt/jarfiles/trimmomatic/v0.39/Trimmomatic-0.39/trimmomatic-0.39.jar'
fastqc_path = '/data/opt/programs/etc/fastqc/FastQC-v0.11.9/fastqc'

# create_directory() accepts one absolute or relative output directory, 'outdir', that makes directory for output files
# if that current directory doesn't exist for the purposes of not accidentally removing a prior directory.
def create_directory(outdir):
    """Function to create the output directory"""
    if os.path.isdir(outdir):
        raise Exception("Directory already exists")
    if not os.path.isdir(outdir):
        os.system("mkdir -p %s" % outdir)
    return

def fastqc_command(perl_path, fastqc_path, threads, sample_name, outdir, read1, read2):
    fastqc_cwd = ['{0}'.format(perl_path), '{0}'.format(fastqc_path), "--threads", threads,
                      "--outdir", '{0}/{1}_fastqc'.format(outdir, sample_name), read1, read2]
    subprocess.run(fastqc_cwd)
    return

# Need to figure out the step parameterization that is best suited for our reads; note that trimmomatic logfile
# stdout is trimming information for each read, not a summary file

def trimmomatic_command(java_path, trimmomatic_path, threads, sample_name, outdir, read1, read2, adapter_database):
    trimmomatic_cwd = ['{0}'.format(java_path), "-jar", '{0}'.format(trimmomatic_path), "PE", "-phred33",
                       "-threads", threads, read1, read2, '{0}/{1}_ILMN_R1.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_unpaired_R1.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_R2.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_unpaired_R2.fastq.gz'.format(outdir, sample_name),
                       'ILLUMINACLIP:{0}:2:30:10'.format(adapter_database), 'LEADING:3', 'TRAILING:15',
                       'SLIDINGWINDOW:4:15', 'MINLEN:36']
    subprocess.run(trimmomatic_cwd)
    return

def kraken2_command(kraken2_path, threads, kraken2_report, sequence, kraken2_database):
    kraken2_cwd = ['{0}'.format(kraken2_path), "--db", kraken2_database, "--confidence", "0.1", "--threads", threads, "--output", '-',
                   "--report", kraken2_report, sequence]
#                   "--output", '{0}/{1}_kraken2_sequence.tsv'.format(outdir, sample_name),
#                   "--report", '{0}/{1}_kraken2_report.tsv'.format(outdir, sample_name), sequence]
    subprocess.run(kraken2_cwd)
    return

def spades_command(spades_path, threads, sample_name, outdir, read1, read2):
    spades_cwd = ('{0}'.format(spades_path), "-t", threads, "--careful", "--cov-cutoff", "auto", "-1", read1,
               "-2", read2, "-o", '{0}/{1}_SPAdes_assembly'.format(outdir, sample_name))
    subprocess.run(spades_cwd)
    return

def contig_filter_command(reformat_path, input, output):
    contig_filter_cwd = ('{0}'.format(reformat_path), 'in={0}'.format(input), 'out={0}'.format(output), 'minlength=500')
    subprocess.run(contig_filter_cwd)
    return

def quast_command(quast_path, threads, sample_name, outdir, input):
    quast_cwd = ('{0}'.format(quast_path), '--threads', threads,
                 '-o', '{0}/{1}_QUAST_results'.format(outdir, sample_name), input)
    subprocess.run(quast_cwd)
    return

def make_prokka_command(prokka_path, contigs, sample_name, outdir, threads):
    prokka_command = ['{0}'.format(prokka_path), '--cpus', threads, '--compliant', '--outdir',
                      '{0}/{1}_prokka_results'.format(outdir, sample_name), '--prefix',
                      '{0}_prokka'.format(sample_name), contigs]
    subprocess.run(prokka_command)
    return

def make_mlst_command(mlst_path, contigs, sample_name, outdir):
    mlst_command = ['{0}'.format(mlst_path), contigs]
    result = subprocess.run(mlst_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    mlst = result.stdout.decode('utf-8')
    mlst_file_name = "{0}/{1}_mlst.tsv".format(outdir, sample_name)
    mlst_file = open(mlst_file_name, 'w')
    mlst_file.write(mlst)
    mlst_file.close()
    return


def make_abricate_command(abricate_path, contigs, db, sample_name, outdir):
    abricate_command = ['{0}'.format(abricate_path), '--nopath', '--mincov', '90', '--minid', '90', '--db', db, contigs]
    result = subprocess.run(abricate_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    abricate = result.stdout.decode('utf-8')
    abricate_file_name = "{0}/{1}_{2}_abricate.tsv".format(outdir, sample_name, db)
    abricate_file = open(abricate_file_name, 'w')
    abricate_file.write(abricate)
    abricate_file.close()
    return

def get_arguments():
    """Parse assembler arguments"""
    parser = argparse.ArgumentParser(description="SPAdes assembly pipeline", add_help=False)

    # Help arguments
    help_group = parser.add_argument_group("Help")
    help_group.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    help_group.add_argument('-V', '--version', action='version', version='%(prog)s version 0.2',
                            help="Show Assembler's version number")

    # input_arguments
    input_group = parser.add_argument_group("Inputs")
    input_group.add_argument('-s', '--sample_name', required=True, help="Name of the sample to use in the "
                                                                        "outdir/outfiles as prefix", type=str,
                             default=None)
    input_group.add_argument('-R1', '--read1', required=True, help='Raw read1 from Illumina sequencing platform', type=str,
                             default=None)
    input_group.add_argument('-R2', '--read2', required=True, help='Raw read2 from Illumina sequencing platform', type=str,
                             default=None)


    # Pipeline arguments
    pipeline_group = parser.add_argument_group("Pipeline Arguments")
    pipeline_group.add_argument('--java_path', required=False, help="Path to java executable. Please use "
                                "\'java\' with path", type=str, default='java')
    pipeline_group.add_argument('--perl_path', required=False, help="Path to perl executable. Please use "
                                "\'perl\' with path", type=str, default='perl')
#    pipeline_group.add_argument('--fastqc_path', required=False, help="Path to java executable. Please use "
#                                "\'fastqc\' with path", type=str, default='fastqc')
#    pipeline_group.add_argument('--trimmomatic_path', required=False, help="Path to trimmomatic executable. Please use "
#                                "\'trimmomatic.jar\' with path", type=str, default='trimmomatic-0.39.jar')
    pipeline_group.add_argument('--kraken2_path', required=False, help="Path to kraken2 executable. Please use "
                                "\'kraken2\' with path", type=str, default='kraken2')
    pipeline_group.add_argument('--spades_path', required=False, help="Path to SPAdes executable. Please use "
                                "\'spades.py\' with path", type=str, default='spades.py')
    pipeline_group.add_argument('--reformat_path', required=False, help="Path to reformat.sh in bbmap suite executable. Please use "
                                "\'reformat.sh\' with path", type=str, default='reformat.sh')
    pipeline_group.add_argument('--quast_path', required=False, help="Path to quast.py script. Please use "
                                "\'quast.py\' with path", type=str, default='quast.py')
    pipeline_group.add_argument('--mlst_path', required=False, help="Path to mlst executable. Please use "
                                "\'mlst\' with path", type=str, default='mlst')
    pipeline_group.add_argument('--prokka_path', required=False, help="Path to prokka executable. Please use "
                                "\'prokka\' with path", type=str,default='prokka')
    pipeline_group.add_argument('--abricate_path', required=False, help="Path to abricate executable. Please use "
                                "\'abricate\' with path", type=str, default='abricate')


    # Optional and Output arguments
    optional_group = parser.add_argument_group("Output and Options")
    optional_group.add_argument('-o', '--outdir', required=True, help="Name of the output directory", type=str,
                                default=None)
    optional_group.add_argument('-t', '--threads', required=False, help="Number of threads to run program", type=str,
                                default='1')

    args = parser.parse_args()
    return args


# Main function that takes argparse arguments and can be passed via a command line prompt.
def run_conditions():
    args = get_arguments()
    outdir1 = args.outdir
    sample_name = args.sample_name
    threads = args.threads
    read1 = args.read1
    read2 = args.read2

    # Load pipeline module
#    module_command(args.module_path, 'spades_pipeline/SPAdes_pipeline-v0.1')
    create_directory(args.outdir)
    start = time.time()
    # perform initial fastqc check on raw reads
    print("\nPerforming quality control checks on input reads\n")
    os.system("mkdir %s" % ('{0}/{1}_fastqc'.format(outdir1, sample_name)))
    fastqc_command(args.perl_path, fastqc_path, threads, sample_name, outdir1, read1, read2)

    # perform trimming using trimmomatic-v0.39
    print("\nTrimming low-quality bases and adapters from raw Illumina short-reads\n")
    trimmomatic_command(args.java_path, trimmomatic_path, threads, sample_name, outdir1, read1, read2, adapter_database1)
    os.system("rm %s" % ('{0}/*unpaired*'.format(outdir1)))

    # perform fastqc check on trimmed reads
    print("\nPerforming quality control checks on trimmed reads\n")
    trim_read1 = '{0}/{1}_ILMN_R1.fastq.gz'.format(outdir1, sample_name)
    trim_read2 = '{0}/{1}_ILMN_R2.fastq.gz'.format(outdir1, sample_name)
    fastqc_command(args.perl_path, fastqc_path, threads, sample_name, outdir1, trim_read1, trim_read2)

    # perform kraken2 contamination check with trimmed read1
    print("\nChecking for contamination in trimmed read1\n")
#    kraken2_trim_R1_output = '{0}/{1}_kraken2_trimR1_sequence.tsv'.format(outdir1, sample_name)
    kraken2_trim_R1_report = '{0}/{1}_kraken2_trimR1_report.tsv'.format(outdir1, sample_name)
    kraken2_command(args.kraken2_path, threads, kraken2_trim_R1_report, trim_read1, kraken2_database1)

    # perform spades-v3.14.0 assembly
    print("\nPerforming SPAdes assembly\n")
    spades_command(args.spades_path, threads, sample_name, outdir1, trim_read1, trim_read2)

    # Can adjust the minLength of contig in the contig_filter_command above
    print("\nFiltering contig lengths by predetermined size\n")
    os.system("cd %s" % ('{0}'.format(outdir1)))
    raw_assembly = './{0}/{1}_SPAdes_assembly/contigs.fasta'.format(outdir1,sample_name)
    filter_assembly = './{0}/{1}_contigs_filtered.fasta'.format(outdir1, sample_name)
    contig_filter_command(args.reformat_path, raw_assembly, filter_assembly)
    os.system("mv %s %s" % ("{0}/{1}_SPAdes_assembly/spades.log".format(outdir1, sample_name),
                            "{0}/{1}_spades.log".format(outdir1, sample_name)))
    os.system('rm -r %s' % '{0}/{1}_SPAdes_assembly'.format(outdir1, sample_name))


    # Perform QUAST quality check on filtered SPAdes assembly
    print("\nPerforming QC on raw, filtered SPAdes assembly\n")
    quast_command(args.quast_path, threads, sample_name, outdir1, filter_assembly)

    # Performing basic preliminary analysis; Need to module load perl-5.28.0 at this moment; probably can fix when
    # Blake can load appropriate modules on 3.28.0 that will not conflict with the perl-5.16 when running fastqc

#    module_command('perl/perl-5.28.0')

    # Perform Prokka
    print("\nCreating annotation files from input fasta file\n")
    make_prokka_command(args.prokka_path, filter_assembly, sample_name, outdir1, threads)
    gbk_file = '{0}/{1}_prokka_results/{1}_prokka.gbk'.format(outdir1, sample_name)
    # Perform MLST, abricate on filtered contigs
    print("\nPerform MLST\n")
    make_mlst_command(args.mlst_path, filter_assembly, sample_name, outdir1)

    print("\nPerform abricate database searches\n")
    make_abricate_command(args.abricate_path, gbk_file, 'card', sample_name, outdir1)
    make_abricate_command(args.abricate_path, gbk_file, 'plasmidfinder', sample_name, outdir1)

    end = time.time()


# Perform some basic checks for assembly
    print("\nTotal Wall Time: ", end - start, "seconds\n")
    print("\nFin! Enjoy your day!\n")

if __name__ == '__main__': run_conditions()
