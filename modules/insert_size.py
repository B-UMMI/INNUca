#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
insert_size.py - Determines the sequencing insert size by reads mapping
<https://github.com/B-UMMI/INNUca/modules/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: October 10, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os.path
from utils import runCommandPopenCommunicate as utils_run_command
from utils import runTime as utils_run_time
from utils import checkPrograms as utils_check_programs
from utils import timer as utils_timer
from functools import partial
import argparse
import time

version = '1.0'


def index_sequence_bowtie2(reference, outdir, threads=1):
    """
    Index reference sequence for Bowtie2

    Parameters
    ----------
    reference : str
        Path to the reference fasta file against which the reads will be mapped
    outdir : str
        Path to the output directory
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if bowtie2-build ran successfully or not
    reference_index : str or None
        Path to the reference Bowtie2 index (if ran successfully, else returns None)
    """

    command = ['bowtie2-build', '--threads', str(threads), reference, os.path.join(outdir, os.path.basename(reference))]
    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=True)

    reference_index = None
    if run_successfully:
        reference_index = os.path.join(outdir, os.path.basename(reference))

    return run_successfully, reference_index


def mapping_bowtie2(fastq, reference_index, outdir, threads=1):
    """
    Map reads against a reference fasta file

    Parameters
    ----------
    fastq : list
        List of fastq files (only two, paired-end reads)
    reference_index : str
        Path to the reference Bowtie2 index
    outdir : str
        Path to the output directory
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Assembly_Mapping module ran successfully or not
    sam : str or None
        If everything went fine, it returns the path to the sam file, otherwise it returns None
    """

    sam = os.path.join(outdir, str('alignment.sam'))

    command = ['bowtie2', '-q', '--very-sensitive', '--threads', str(threads), '-x', reference_index, '-1', fastq[0],
               '-2', fastq[1], '--no-mixed', '--no-discordant', '--no-unal', '-S', sam]

    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=True)

    if not run_successfully:
        sam = None

    return run_successfully, sam


def get_statistics_samtools(alignment, outdir):
    """
    Run Samtools stats to get several statistics from the alignment file

    Parameters
    ----------
    alignment : str
        Path to the alignment file (can be SAM, BAM or CRAM)
    outdir : str
        Path to the output directory

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Assembly_Mapping module ran successfully or not
    samtools_stats : str or None
        If everything went fine, it returns the path to the samtools stats file, otherwise it returns None
    """

    samtools_stats = os.path.join(outdir, 'samtools_stats.txt')

    command = ['samtools', 'stats', alignment, '>', samtools_stats]

    run_successfully, _, _ = utils_run_command(command=command, shell_True=True, timeout_sec_None=None,
                                               print_comand_True=True)

    if not run_successfully:
        samtools_stats = None

    return run_successfully, samtools_stats


def parse_statistics_samtools(samtools_stats):
    """
    Parse Samtools statistics and get useful information (insert_size, insert_size_sd)

    Parameters
    ----------
    samtools_stats : str
        Path to the samtools stats file

    Returns
    -------
    statistics : dict
        Dictionary with the statistics. Keys are statistics description and values are the correspondent values
    """

    statistics = {}
    with open(samtools_stats, 'rt') as reader:
        counter = 0
        for line in reader:
            if counter <= 2:
                if line.startswith('SN'):
                    line = line.split('\t')
                    if line[1] == 'insert size average:':
                        statistics['insert_size'] = round(float(line[2]), 0)
                        counter += 1
                    elif line[1] == 'insert size standard deviation:':
                        statistics['insert_size_sd'] = float(line[2])
                        counter += 1
            else:
                break

    return statistics


def write_reports(statistics, outdir):
    """
    Write Samtools stats report

    Parameters
    ----------
    statistics : dict
        Dictionary with the statistics. Keys are statistics description and values are the correspondent values
    outdir : str
        Path to the output directory

    Returns
    -------

    """

    with open(os.path.join(outdir, 'samtools_stats_report.tab'), 'wt') as writer:
        writer.write('#' + '\t'.join(sorted(statistics.keys())) + '\n')
        writer.write('\t'.join([str(statistics[k]) for k in sorted(statistics.keys())]) + '\n')


def clean_outdir(outdir, reference_index=None, alignment=None):
    """
    Remove temporary files

    Parameters
    ----------
    outdir : str
        Path to the output directory
    reference_index : str or None
        Path to the reference Bowtie2 index (or None if only using alignment file for the analysis)
    alignment : str or None
        Path to the alignment file (or None if don't want to remove the alignment file)

    Returns
    -------

    """

    files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
    for file_found in files:
        # Bowtie2 index
        if reference_index is not None:
            if file_found.startswith(os.path.basename(reference_index) + '.') and file_found.endswith('.bt2'):
                os.remove(os.path.join(outdir, file_found))
        # Alignment
        if alignment is not None:
            if file_found == os.path.basename(reference_index):
                os.remove(os.path.join(outdir, file_found))


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 insert_size.py"')

    parser = argparse.ArgumentParser(prog='insert_size.py',
                                     description='Determines the sequencing insert size by reads mapping',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-r', '--reference', nargs=1, type=argparse.FileType('r'), required=True,
                                 metavar='/path/to/reference.fasta',
                                 help='Path to the reference fasta file. Ideally it must be the sample assembly.')
    parser_required.add_argument('-f', '--fastq', nargs=2, type=argparse.FileType('r'), required=True,
                                 metavar=('/path/to/fastq/file_1.fq.gz', '/path/to/fastq/file_2.fq.gz'),
                                 help='Path to the paired-end fastq files')

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the sequences will be stored (default: ./)',
                                         required=False, default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N',
                                         help='Number of threads to use (default: 1)', required=False, default=1)

    args = parser.parse_args()

    start_time = time.time()

    print('\n'
          '===>  RUNNING  insert_size.py  <===\n')

    missing_programs, _ = utils_check_programs({'bowtie2': ['--version', '>=', '2.2.9'],
                                                'samtools': ['--version', '==', '1.3.1']})
    if len(missing_programs) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missing_programs))

    print()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    args.reference = os.path.abspath(args.reference[0].name)
    args.fastq = [os.path.abspath(fastq.name) for fastq in args.fastq]

    run_successfully, reference_index = index_sequence_bowtie2(reference=args.reference, outdir=args.outdir,
                                                               threads=args.threads)

    sam = None
    if run_successfully:
        run_successfully, sam = mapping_bowtie2(fastq=args.fastq, reference_index=reference_index, outdir=args.outdir,
                                                threads=args.threads)

        if run_successfully:
            run_successfully, samtools_stats = get_statistics_samtools(alignment=sam, outdir=args.outdir)

            if run_successfully:
                statistics = parse_statistics_samtools(samtools_stats=samtools_stats)
                write_reports(statistics=statistics, outdir=args.outdir)

    clean_outdir(outdir=args.outdir, reference_index=reference_index, alignment=sam)

    _ = utils_run_time(start_time)

    if not run_successfully:
        sys.exit('Something went wrong!')


if __name__ == "__main__":
    main()
