#!/usr/bin/env python2

# -*- coding: utf-8 -*-

"""
trueCoverage_rematch.py - Estimate the true bacterial chromosome
coverage and detects contamination with different strain or species.

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: February 27, 2018

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

import os.path
import utils
import functools
import sys
import argparse

version = '0.2'


def check_existing_default_config(species, script_path):
    species = species.lower().split(' ')
    true_coverage_config_folder = os.path.join(os.path.dirname(script_path), 'modules', 'trueCoverage_rematch', '')
    config = None
    reference = None
    files = [f for f in os.listdir(true_coverage_config_folder) if not f.startswith('.') and
             os.path.isfile(os.path.join(true_coverage_config_folder, f))]
    for file_found in files:
        file_path = os.path.join(true_coverage_config_folder, file_found)
        if file_found == '_'.join(species) + '.config':
            config = file_path
        elif file_found == '_'.join(species) + '.fasta':
            reference = file_path
    return config, reference


def parse_config(config_file):
    config = {'reference_file': None,
              'length_extra_seq': None,
              'maximum_number_absent_genes': None,
              'maximum_number_genes_multiple_alleles': None,
              'minimum_read_coverage': None,
              'minimum_depth_presence': None,
              'minimum_depth_call': None,
              'minimum_depth_frequency_dominant_allele': None,
              'minimum_gene_coverage': None,
              'minimum_gene_identity': None}

    with open(config_file, 'rtU') as reader:
        field = None
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                line = line.split(' ')[0]
                if line.startswith('#'):
                    line = line[1:].split(' ')[0]
                    field = line
                else:
                    if field is not None:
                        if field in ['length_extra_seq', 'maximum_number_absent_genes',
                                     'maximum_number_genes_multiple_alleles', 'minimum_read_coverage',
                                     'minimum_depth_presence', 'minimum_depth_call', 'minimum_gene_coverage',
                                     'minimum_gene_identity']:
                            line = int(line)
                            if field in ['minimum_gene_coverage', 'minimum_gene_identity']:
                                if line < 0 or line > 100:
                                    sys.exit('minimum_gene_coverage in trueCoverage_rematch config file must be an'
                                             ' integer between 0 and 100')
                        elif field == 'minimum_depth_frequency_dominant_allele':
                            line = float(line)
                            if line < 0 or line > 1:
                                sys.exit('minimum_depth_frequency_dominant_allele in trueCoverage_rematch config file'
                                         ' must be a double between 0 and 1')
                        config[field] = line
                        field = None

    for field in config:
        if config[field] is None:
            sys.exit(field + ' in trueCoverage_rematch config file is missing')

    return config


def clean_headers_reference_file(reference_file, outdir, extra_seq, rematch_module):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    print 'Checking if reference sequences contain ' + str(problematic_characters) + '\n'
    new_reference_file = str(reference_file)
    sequences, genes, headers_changed = rematch_module.get_sequence_information(reference_file, extra_seq)
    if headers_changed:
        print 'At least one of the those characters was found. Replacing those with _' + '\n'
        new_reference_file = os.path.join(outdir,
                                          os.path.splitext(os.path.basename(reference_file))[0] +
                                          '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = rematch_module.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, genes, sequences


def rematch_report_assess_failing(outdir, time_str, rematch_folder, sample_data_general, config):
    """
    Copy ReMatCh sample report to outdir and assess sample failing status

    Parameters
    ----------
    outdir : str
        Path to the directory where true_coverage results will be stored, e.g. "/path/to/output/directory/"
    time_str : str or None
        String containing date and time information in the following format "%Y%m%d-%H%M%S" or None
    rematch_folder : str
        Path to the temporary ReMatCh folder that contain ReMatCh results
    sample_data_general : dict
        ReMatCh sample_data_general dictionary containing general sample results
    config : dict
        parse_config config dictionary containing the settings to run trueCoverage_rematch

    Returns
    -------
    failing : dict
        Dictionary containing the reasons for failing true_coverage
    """

    print 'Writing report file'
    os.rename(os.path.join(rematch_folder, 'rematchModule_report.txt'),
              os.path.join(outdir, 'trueCoverage_report.{time_str}.txt'.format(time_str=time_str)
              if time_str is not None else 'trueCoverage_report.txt'))

    failing = {}
    if sample_data_general['number_absent_genes'] > config['maximum_number_absent_genes']:
        failing['absent_genes'] = 'The number of absent genes ({real_absent}) exceeds the maximum allowed' \
                                  ' ({max_absent})'.format(real_absent=sample_data_general['number_absent_genes'],
                                                           max_absent=config['maximum_number_absent_genes'])
    if sample_data_general['number_genes_multiple_alleles'] > config['maximum_number_genes_multiple_alleles']:
        failing['multiple_alleles'] = 'The number of genes with multiple alleles' \
                                      ' ({real_multiple}) exceeds the maximum allowed' \
                                      ' ({max_multiple})'.format(
            real_multiple=sample_data_general['number_genes_multiple_alleles'],
            max_multiple=config['maximum_number_genes_multiple_alleles'])
    if sample_data_general['mean_sample_coverage'] < config['minimum_read_coverage']:
        failing['read_coverage'] = 'The mean read coverage for genes present' \
                                   ' ({real_coverage}) dit not meet the minimum required' \
                                   ' ({min_coverage})'.format(real_coverage=sample_data_general['mean_sample_coverage'],
                                                              min_coverage=config['minimum_read_coverage'])
    return failing


trueCoverage_timer = functools.partial(utils.timer, name='True coverage check')


@trueCoverage_timer
def runTrueCoverage(sample, fastq, reference, threads, outdir, extra_seq, min_cov_presence, min_cov_call,
                    min_frequency_dominant_allele, min_gene_coverage, debug, min_gene_identity,
                    true_coverage_config, rematch_script, conserved_true=True, num_map_loc=1):
    pass_qc = False
    failing = {}

    true_coverage_folder = os.path.join(outdir, 'trueCoverage', '')
    utils.removeDirectory(true_coverage_folder)
    os.mkdir(true_coverage_folder)

    sys.path.append(os.path.join(os.path.dirname(rematch_script), 'modules', ''))
    import rematch_module

    # Run ReMatCh
    reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(reference, true_coverage_folder,
                                                                                       extra_seq, rematch_module)
    time_taken, run_successfully, data_by_gene, sample_data_general, consensus_files, consensus_sequences = rematch_module.runRematchModule(sample, fastq, reference_file, threads, true_coverage_folder, extra_seq, min_cov_presence, min_cov_call, min_frequency_dominant_allele, min_gene_coverage, conserved_true, debug, num_map_loc, min_gene_identity, 'first', 7, 'none', reference_dict, 'X', None, gene_list_reference, True)

    if run_successfully:
        failing = rematch_report_assess_failing(outdir, None, true_coverage_folder, sample_data_general, true_coverage_config)
    else:
        failing['sample'] = 'Did not run'

    if len(failing) == 0:
        pass_qc = True
        failing['sample'] = False
    else:
        print failing

    if not debug:
        utils.removeDirectory(true_coverage_folder)

    return run_successfully, pass_qc, failing


def arguments_required_length(tuple_length_options, argument_name):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) not in tuple_length_options:
                msg = 'argument {argument_name} requires one of the following number of arguments:' \
                      ' {tuple_length_options}'.format(argument_name=argument_name,
                                                       tuple_length_options=tuple_length_options)
                parser.error(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def arguments_check_directory(argument_name):
    """
    Check if the directory passed to the argument exists

    Parameters
    ----------
    argument_name : str
        Argument name, e.g. '--indir'

    Returns
    -------
    ArgumentsCheckDirectory : str
        Full path of the directory
    """
    class ArgumentsCheckDirectory(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            directory = os.path.abspath(values)
            if not os.path.isdir(directory):
                msg = 'argument {argument_name}: {directory} is not a directory'.format(argument_name=argument_name,
                                                                                        directory=directory)
                parser.error(msg)
            setattr(args, self.dest, directory)
    return ArgumentsCheckDirectory


def check_fasta_config_exist(indir, species):
    """
    Check if species fasta and config files exist inside indir

    Parameters
    ----------
    indir : str
        Directory path
    species : list
        List with species name, e.g. ['escherichia', 'coli']

    Returns
    -------
    msg : str or None
        Message with the error found or None
    required_files : dict
        Dictionary with the two files required, e.g.
        {'fasta': /path/escherichia_coli.fasta, 'config': /path/escherichia_coli.config}
    """

    files = [f for f in os.listdir(indir) if not f.startswith('.') and os.path.isfile(os.path.join(indir, f))]
    required_files = {}
    for file_found in files:
        root, ext = os.path.splitext(file_found)
        if root.lower() == '_'.join(species).lower():
            ext = ext.lstrip('.')
            if ext in ('fasta', 'config'):
                required_files[ext] = os.path.join(indir, file_found)
    msg = None
    if len(required_files) == 1:
        msg = 'only found the {ext} file for {species} species (both fasta and config are' \
              ' required)'.format(ext=required_files.keys()[0], species=' '.join(species))
    elif len(required_files) == 0:
        msg = 'no files were found for {species} species (a fasta and config files are' \
              ' required)'.format(species=' '.join(species))

    return msg, required_files


def import_rematch(minimum_rematch_version):
    """
    Check if ReMatCh is in the PATH and has the correct version

    Parameters
    ----------
    minimum_rematch_version : str
        Minimum ReMatCh version, e.g. "3.2"

    Returns
    -------
    rematch_module : python module
        If everything is OK, it returns rematch_module python module
    """

    import subprocess

    command = ['which', 'rematch.py']
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = p.communicate()
    if p.returncode == 0:
        rematch_script = stdout.rstrip('\r\n')
        command = ['rematch.py', '--version']
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, stderr = p.communicate()
        if p.returncode == 0:
            stderr = stderr.rstrip('\r\n').split(' ')[-1].replace('v', '')
            print('rematch.py ({version}) found'.format(version=stderr))
            program_found_version = stderr.split('.')
            minimum_rematch_version = minimum_rematch_version.split('.')
            if len(minimum_rematch_version) == 3:
                if len(program_found_version) == 2:
                    program_found_version.append(0)
                else:
                    program_found_version[2] = program_found_version[2].split('_')[0]
            for i in range(0, len(minimum_rematch_version)):
                if int(program_found_version[i]) > int(minimum_rematch_version[i]):
                    break
                elif int(program_found_version[i]) == int(minimum_rematch_version[i]):
                    continue
                else:
                    sys.exit('It is required at least ReMatCh with version'
                             ' v{minimum_rematch_version}'.format(minimum_rematch_version=minimum_rematch_version))

            sys.path.append(os.path.join(os.path.dirname(rematch_script), 'modules', ''))
            import rematch_module

            bowtie2 = os.path.join(os.path.dirname(rematch_script), 'src', 'bowtie2-2.2.9')
            samtools = os.path.join(os.path.dirname(rematch_script), 'src', 'samtools-1.3.1', 'bin')
            bcftools = os.path.join(os.path.dirname(rematch_script), 'src', 'bcftools-1.3.1', 'bin')
            os.environ['PATH'] = str(':'.join([bowtie2, samtools, bcftools, os.environ['PATH']]))
            print('PATH={path}'.format(path=os.environ['PATH']))

            return rematch_module
        else:
            sys.exit('It was not possible to determine ReMatCh version')
    else:
        sys.exit('ReMatCh not found in PATH')


def start_logger(workdir):
    import time
    time_str = time.strftime("%Y%m%d-%H%M%S")
    sys.stdout = Logger(workdir, time_str)
    logfile = sys.stdout.getLogFile()
    return logfile, time_str


class Logger(object):
    def __init__(self, out_directory, time_str):
        self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
        self.terminal = sys.stdout
        self.log = open(self.logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        pass

    def getLogFile(self):
        return self.logfile


def main():
    if sys.version_info[0] > 2:
        sys.exit('Must be using Python 2. Try calling "python2 trueCoverage_rematch.py"')

    parser = argparse.ArgumentParser(prog='python2 trueCoverage_rematch.py',
                                     description="Estimate the true bacterial chromosome coverage and detects"
                                                 " contamination with different strain or species",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-f', '--fastq', nargs='+', action=arguments_required_length((1, 2), '--fastq'),
                                 type=argparse.FileType('r'), metavar=('/path/to/input/file.fq.gz'),
                                 help='Path to single OR paired-end fastq files. If two files are passed, they will be'
                                      ' assumed as being the paired fastq files',
                                 required=True)

    parser_true_coverage_files = parser.add_argument_group('Options for trueCoverage_rematch species files')
    parser_true_coverage_files.add_argument('-s', '--species', nargs=2, type=str,
                                            metavar=('Yersinia', 'enterocolitica'), help='Species name',
                                            required=False)
    parser_true_coverage_files.add_argument('-i', '--indir', type=str, action=arguments_check_directory('--indir'),
                                            metavar='/path/to/fasta/config/indir/directory/',
                                            help='Path to the directory where species reference fasta files and config'
                                                  ' files can be found',
                                            required=False)

    parser_external_files = parser.add_argument_group('Options for external fasta and config files')
    parser_external_files.add_argument('-r', '--reference', type=argparse.FileType('r'),
                                       metavar='/path/to/reference_sequence.fasta',
                                       help='Fasta file containing reference sequences. Ideally they should be'
                                            ' housekeeping gene sequences distributed throughout the genome.'
                                            ' Alternativelly, MLST gene sequences might be a good approximation.',
                                       required=False)
    parser_external_files.add_argument('-c', '--config', type=argparse.FileType('r'), metavar='/path/to/config.file',
                                       help='Config file with the settings to run trueCoverage_rematch. Check some'
                                            ' examples in INNUca GitHub ('
                                            'https://github.com/B-UMMI/INNUca/tree/master/modules/trueCoverage_rematch).',
                                       required=False)

    parser_optional = parser.add_argument_group('Facultative options')
    parser_optional.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                 help='Path to the directory where the results will be stored',
                                 required=False, default='.')
    parser_optional.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use',
                                 required=False, default=1)
    parser_optional.add_argument('--json', action='store_true', help='Stores the results also in JSON format')

    args = parser.parse_args()

    if (args.species is not None or args.indir is not None) and (args.reference is not None or args.config is not None):
        parser.error('Do not mix options from the two pairs: --species and --indir OR --reference and --config')
    elif args.species is None and args.indir is None and args.reference is None and args.config is None:
        parser.error('At least one of the following option pairs must be specified:'
                     ' --species and --indir OR --reference and --config')

    if args.species is not None:
        error_msg, required_files = check_fasta_config_exist(args.indir, args.species)
        if error_msg is not None:
            parser.error('argument {argument_name}: {error_msg}'.format(argument_name='--indir', error_msg=error_msg))
    else:
        required_files = {'fasta': os.path.abspath(args.reference.name), 'config': os.path.abspath(args.config.name)}

    rematch_module = import_rematch('3.2')

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # Start logger
    logfile, time_str = start_logger(args.outdir)

    config = parse_config(required_files['config'])

    import tempfile
    rematch_folder = tempfile.mkdtemp(prefix='trueCoverage_rematch_', suffix='_tmp', dir=args.outdir)

    reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(required_files['fasta'],
                                                                                       rematch_folder,
                                                                                       config['length_extra_seq'],
                                                                                       rematch_module)
    time_taken, run_successfully, data_by_gene, sample_data_general, consensus_files, consensus_sequences = rematch_module.runRematchModule('sample', [os.path.abspath(fastq.name) for fastq in args.fastq], reference_file, args.threads, rematch_folder, config['length_extra_seq'], config['minimum_depth_presence'], config['minimum_depth_call'], config['minimum_depth_frequency_dominant_allele'], config['minimum_gene_coverage'], True, True, 1, config['minimum_gene_identity'], 'first', 7, 'none', reference_dict, 'X', None, gene_list_reference, True)

    import shutil
    if run_successfully:
        failing = rematch_report_assess_failing(args.outdir, time_str, rematch_folder, sample_data_general, config)
        if len(failing) > 0:
            with open(os.path.join(args.outdir, 'failing.' + time_str + '.txt'), 'wt') as writer:
                for scope, reason in failing.items():
                    writer.write('#{scope}\n'
                                 '{reason}\n'.format(scope=scope, reason=reason))
        if args.json:
            import json
            with open(os.path.join(args.outdir, 'sample_data_general.' + time_str + '.json'), 'wt') as writer:
                json.dump(sample_data_general, writer, separators=(",", ":"))
            if len(failing) > 0:
                with open(os.path.join(args.outdir, 'failing.' + time_str + '.json'), 'wt') as writer:
                    json.dump(failing, writer, separators=(",", ":"))
    else:
        shutil.rmtree(rematch_folder)
        sys.exit('Something went wrong while running ReMatCh')

    shutil.rmtree(rematch_folder)


if __name__ == "__main__":
    main()
