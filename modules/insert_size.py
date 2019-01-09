#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
insert_size.py - Determines the sequencing insert size by reads mapping
<https://github.com/B-UMMI/INNUca/modules/>

Copyright (C) 2019 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: January 02, 2019

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
from utils import removeDirectory as utils_remove_directory
from functools import partial
import argparse
import time

try:
    import plotly.offline as plot_off
    import plotly.graph_objs as graph_obj
except ImportError as e:
    print('WARNING: {}'.format(e))

version = '1.3'


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
    print('')

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

    command = ['bowtie2', '-q', '--very-fast', '--threads', str(threads), '-x', reference_index, '-1', fastq[0],
               '-2', fastq[1], '--fr', '-I', '0', '-X', '2000', '--no-discordant', '--no-mixed', '--no-unal', '-S', sam]

    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=True)
    print('')

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
    print('')

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

    with open(os.path.join(outdir, 'insert_size_report.tab'), 'wt') as writer:
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
            if file_found == os.path.basename(alignment):
                os.remove(os.path.join(outdir, file_found))


def prepare_insert_size_distribution(samtools_stats):
    """
    Collect the data to produce a distribution plot of the insert sizes

    Parameters
    ----------
    samtools_stats : str
        Path to the samtools stats file

    Returns
    -------
    x : list
        List of x axis values
    y : list
        List of y axis values
    """

    x = []
    y = []
    with open(samtools_stats, 'rt') as reader:
        for line in reader:
            if line.startswith('IS'):
                line = line.split('\t')
                if int(line[2]) > 0:
                    x.append(int(line[1]))
                    y.append(int(line[2]))

    return x, y


def absolute_counts_2_frequencies(counts):
    """
    Converts absolute counts to frequencies

    Parameters
    ----------
    counts : list
        List of absolute counts that will be converted to frequencies

    Returns
    -------
    normalized_values : list
        List of absolute counts already converted to frequencies
    """
    normalized_values = list(map(lambda x: float(x) / sum(counts), counts))
    return normalized_values


def create_plot_trace_insert_size_distribution(trace_name, x_values, y_values, color='rgb(0, 0, 0)'):
    """
    Collect the data to produce a distribution plot of the insert sizes

    Parameters
    ----------
    trace_name : str
        Name of the trace
    x_values : list
        List of x axis values
    y_values : list
        List of y axis values
    color : str, default 'rgb(0, 0, 0)'
        String with the colour to be used. Something like 'rgb(205, 12, 24)'. Default is black

    Returns
    -------
    trace : graph_obj.Scatter
        plotly.graph_objs.Scatter
        A plot trace for the given dataset
    """

    trace = graph_obj.Scatter(name=trace_name,
                              x=x_values,
                              y=y_values,
                              mode='lines',
                              line=dict(color=color
                                        )
                              )

    return trace


def plot_insert_size_distribution(traces_list, outdir):
    """
    Collect the data to produce a distribution plot of the insert sizes

    Parameters
    ----------
    traces_list : list
        List of plot traces (plotly.graph_objs.Scatter)
    outdir : str
        Path to the output directory

    Returns
    -------
    """

    plot_off.plot({"data": traces_list,
                   "layout": graph_obj.Layout(title="Insert size distribution",
                                              xaxis=dict(title='Insert size'),
                                              yaxis=dict(title='Frequency'))
                   },
                  show_link=False,
                  output_type='file',
                  filename=os.path.join(outdir, 'insert_size_distribution.html'),
                  include_plotlyjs=True,
                  auto_open=False,
                  )


insert_timer = partial(utils_timer, name='Insert size determination')


@insert_timer
def run_for_innuca(sample_name, reference, fastq, outdir, threads=1, dist=False):
    """
    Runs insert_size for INNUca

    Parameters
    ----------
    sample_name : str
        Sample's name
    reference : str
        Path to the sample assembly fasta file
    fastq : list
        List with the paths to the paired-end fastq files
    outdir : str
        Path to the output directory
    threads : int, default 1
        Number of threads to be used
    dist : bool, default False
        True if want to produce a distribution plot of the insert sizes (requires Plotly), else False

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Kraken module ran successfully or not
    pass_qc : None
        QA/QC not performed
    time_taken : float
        Seconds that run_for_innuca took to run
    failing : dict
        Dictionary with the failing reasons. If sample did not fail, it is only {'sample': False}. If it failed, keys
        will be the level of failing, and values list of strings
    """

    failing = {'sample': False}

    module_folder = os.path.join(outdir, 'insert_size', '')
    utils_remove_directory(module_folder)
    os.mkdir(module_folder)

    run_successfully, reference_index = index_sequence_bowtie2(reference=reference, outdir=module_folder,
                                                               threads=threads)

    alignment_file = None
    if run_successfully:
        run_successfully, alignment_file = mapping_bowtie2(fastq=fastq, reference_index=reference_index,
                                                           outdir=module_folder, threads=threads)
        if run_successfully:
            run_successfully, samtools_stats = get_statistics_samtools(alignment=alignment_file, outdir=module_folder)
            if run_successfully:
                statistics = parse_statistics_samtools(samtools_stats=samtools_stats)
                write_reports(statistics=statistics, outdir=outdir)
                if dist:
                    x_values, y_values = prepare_insert_size_distribution(samtools_stats=samtools_stats)
                    y_values = absolute_counts_2_frequencies(y_values)
                    plot_trace = create_plot_trace_insert_size_distribution(trace_name=sample_name, x_values=x_values,
                                                                            y_values=y_values, color='rgb(0, 0, 0)')
                    plot_insert_size_distribution(traces_list=[plot_trace], outdir=outdir)

    clean_outdir(outdir=module_folder, reference_index=reference_index, alignment=alignment_file)

    if not run_successfully:
        failing['sample'] = ['Did not run']
        print(failing['sample'])

    return run_successfully, None, failing


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

    parser_output_options = parser.add_argument_group('Output options')
    parser_output_options.add_argument('--dist', action='store_true',
                                       help='Produces a distribution plot of the insert sizes (requires Plotly)')

    args = parser.parse_args()

    msg = []
    if args.insertSizeDist:
        try:
            import plotly
        except ImportError as import_error:
            msg.append('Used --insertSizeDist but Python Plotly package was not found: {}'.format(import_error))

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

    start_time = time.time()

    print('\n'
          '===>  RUNNING  insert_size.py  <===')

    missing_programs, _ = utils_check_programs({'bowtie2': ['--version', '>=', '2.2.9'],
                                                'samtools': ['--version', '==', '1.3.1']})
    if len(missing_programs) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missing_programs))

    print('')

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
                if args.dist:
                    x_values, y_values = prepare_insert_size_distribution(samtools_stats=samtools_stats)
                    y_values = absolute_counts_2_frequencies(y_values)
                    plot_trace = create_plot_trace_insert_size_distribution(trace_name='', x_values=x_values,
                                                                            y_values=y_values, color='rgb(0, 0, 0)')
                    plot_insert_size_distribution(traces_list=[plot_trace], outdir=args.outdir)

    clean_outdir(outdir=args.outdir, reference_index=reference_index, alignment=sam)

    print('')

    _ = utils_run_time(start_time)

    if not run_successfully:
        sys.exit('Something went wrong!')


if __name__ == "__main__":
    main()
