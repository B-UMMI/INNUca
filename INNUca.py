#!/usr/bin/env python2

# -*- coding: utf-8 -*-


"""
INNUca - Reads Control and Assembly
INNUca.py - INNUENDO quality control of reads, de novo assembly and contigs
quality assessment, and possible contamination detection
<https://github.com/B-UMMI/INNUca>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: February 05, 2020

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

import modules.utils as utils
from modules.kraken import run_for_innuca as kraken
from modules.kraken import get_kraken_version as kraken_version
import modules.fastQintegrity as fastQintegrity
import modules.estimated_coverage as coverage
import modules.fastqc as fastqc
import modules.trimmomatic as trimmomatic
import modules.pear as pear
import modules.spades as spades
import modules.assembly_mapping as assembly_mapping
import modules.pilon as pilon
import modules.mlst as mlst
from modules.insert_size import run_for_innuca as insert_size
import modules.trueCoverage_rematch as trueCoverage
import modules.combine_reports as combine_reports
import time
import os
import sys


# TODO: parse breseq: https://stackoverflow.com/questions/2870667/how-to-convert-an-html-table-to-an-array-in-python


def get_trueCoverage_config(skipTrueCoverage, trueConfigFile, speciesExpected, script_path):
    trueCoverage_config = None
    if not skipTrueCoverage:
        trueCoverage_reference = None

        if trueConfigFile is None:
            print('No trueCoverage_ReMatCh config file was provided. Search for default files')
            trueCoverage_config_file, trueCoverage_reference = trueCoverage.check_existing_default_config(speciesExpected, script_path)
        else:
            trueCoverage_config_file = trueConfigFile

        if trueCoverage_config_file is not None:
            trueCoverage_config = trueCoverage.parse_config(trueCoverage_config_file)
        if trueConfigFile is None and trueCoverage_config is not None:
            trueCoverage_config['reference_file'] = trueCoverage_reference

        if trueCoverage_config is not None:
            print('The following trueCoverage_ReMatCh config file will be used: ' + trueCoverage_config_file)
            print('The following trueCoverage_ReMatCh reference file will be'
                  ' used: {reference}\n'.format(reference=trueCoverage_config['reference_file']))
        else:
            print('No trueCoverage_ReMatCh config file was found')
    return trueCoverage_config


def include_rematch_dependencies_path(do_not_use_provided_software):
    rematch_script = utils.find_rematch()

    if not do_not_use_provided_software and rematch_script is not None:
        path_variable = os.environ['PATH']
        script_folder = os.path.dirname(rematch_script)
        bcftools = os.path.join(script_folder, 'src', 'bcftools-1.3.1', 'bin')
        os.environ['PATH'] = str(':'.join([bcftools, path_variable]))

    return rematch_script


version_kraken_global = None


def main():
    version = '4.2.2'
    args = utils.parseArguments(version)

    general_start_time = time.time()
    time_str = time.strftime("%Y%m%d-%H%M%S")

    # Check if output directory exists
    outdir = os.path.abspath(os.path.join(args.outdir, ''))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Start logger
    if not args.noLog:
        sys.stdout = utils.Logger(outdir, time_str)

    print('\n' + '==========> INNUca.py <==========')
    print('\n' + 'Program start: ' + time.ctime())

    # Tells where the logfile will be stored
    if not args.noLog:
        print('\n' + 'LOGFILE:')
        print(sys.stdout.getLogFile())

    # Print command
    print('\n' + 'COMMAND:')
    script_path = os.path.abspath(sys.argv[0])
    print(sys.executable + ' ' + script_path + ' ' + ' '.join(sys.argv[1:]))

    # Print directory where programme was lunch
    print('\n' + 'PRESENT DIRECTORY:')
    print(os.getcwd())

    # Print program version
    print('\n' + 'VERSION INNUca.py:')
    utils.script_version_git(version=version, current_directory=os.getcwd(), script_path=script_path,
                             no_git_info=args.noGitInfo)

    # Get CPU information
    utils.get_cpu_information(outdir, time_str)

    # Get trueCoverage_ReMatCh settings
    trueCoverage_config = get_trueCoverage_config(args.skipTrueCoverage,
                                                  args.trueConfigFile.name if args.trueConfigFile is not None else None,
                                                  args.speciesExpected, script_path)

    # Check programms
    programs_version_dictionary = {}
    programs_version_dictionary['gunzip'] = {'required': ['--version', '>=', '1.6']}

    # Java check first for java dependents check next
    if not (args.skipFastQC and args.skipTrimmomatic and (args.skipPilon or args.skipSPAdes)):
        # programs_version_dictionary['java'] = ['-version', '>=', '1.8']
        programs_version_dictionary['java'] = {'required': [None, '>=', '1.8']}  # For OpenJDK compatibility
    missingPrograms, programs_version_dictionary = utils.checkPrograms(programs_version_dictionary)
    if len(missingPrograms) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))

    rematch_script = None

    if args.runKraken:
        global version_kraken_global
        version_kraken_global = kraken_version()
        if version_kraken_global == 2:
            programs_version_dictionary['kraken2'] = {'required': ['--version', '>=', '2.0.6']}
        else:
            programs_version_dictionary['kraken'] = {'required': ['--version', '>=', '0.10.6']}
            programs_version_dictionary['kraken-repor'] = {'required': ['--version', '>=', '0.10.6']}
    if not args.skipTrueCoverage and trueCoverage_config is not None:
        rematch_script = include_rematch_dependencies_path(args.doNotUseProvidedSoftware)
        programs_version_dictionary['rematch.py'] = {'required': ['--version', '>=', '4.0.1']}
        programs_version_dictionary['bcftools'] = {'required': ['--version', '==', '1.3.1']}
    if not (args.skipTrueCoverage and ((args.skipAssemblyMapping and args.skipPilon) or args.skipSPAdes)):
        programs_version_dictionary['bowtie2'] = {'required': ['--version', '>=', '2.2.9']}
        programs_version_dictionary['samtools'] = {'required': ['--version', '==', '1.3.1']}
    if not args.skipFastQC:
        programs_version_dictionary['fastqc'] = {'required': ['--version', '==', '0.11.5']}
    if not args.skipTrimmomatic:
        programs_version_dictionary['trimmomatic-{version}.jar'.format(version=args.trimVersion)] = \
            {'required': ['-version', '==', args.trimVersion]}
    if args.runPear:
        programs_version_dictionary['pear'] = {'required': ['--version', '>=', '0.9.10']}
    if not args.skipSPAdes:
        programs_version_dictionary['spades.py'] = {'required': ['--version', '>=', '3.9.0']}
    if not (args.skipPilon or args.skipSPAdes):
        programs_version_dictionary['pilon-{version}.jar'.format(version=args.pilonVersion)] = \
            {'required': ['--version', '==', args.pilonVersion]}
    if not (args.skipMLST or args.skipSPAdes):
        programs_version_dictionary['mlst'] = {'required': ['--version', '>=', '2.4']}
    if args.runInsertSize and not args.skipSPAdes:
        if args.skipAssemblyMapping and args.skipPilon:
            programs_version_dictionary['bowtie2'] = {'required': ['--version', '>=', '2.2.9']}
            programs_version_dictionary['samtools'] = {'required': ['--version', '==', '1.3.1']}

    # Set and print PATH variable
    utils.setPATHvariable(args, script_path)

    missingPrograms, programs_version_dictionary = utils.checkPrograms(programs_version_dictionary)
    if len(missingPrograms) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))

    # .jar paths
    jar_path_trimmomatic = None
    if not args.skipTrimmomatic:
        jar_path_trimmomatic = \
            programs_version_dictionary['trimmomatic-{version}.jar'.format(version=args.trimVersion)]['found']['path']

    jar_path_pilon = None
    if not args.skipPilon and not args.skipSPAdes:
        jar_path_pilon = \
            programs_version_dictionary['pilon-{version}.jar'.format(version=args.pilonVersion)]['found']['path']

    # Get SPAdes version
    spades_version = None
    if not args.skipSPAdes:
        spades_version = programs_version_dictionary['spades.py']['found']['version']

    # pairEnd_filesSeparation_list = args.pairEnd_filesSeparation
    pairEnd_filesSeparation_list = None
    samples, inputDirectory, removeCreatedSamplesDirectories, indir_same_outdir = \
        get_samples(args.inputDirectory, args.fastq, outdir, pairEnd_filesSeparation_list)

    # Start running the analysis
    print('\n' + 'RUNNING INNUca.py')

    # Prepare run report file
    samples_report_path = os.path.join(outdir, 'samples_report.' + time_str + '.tab')
    utils.start_sample_report_file(samples_report_path)

    number_samples_successfully = 0
    number_samples_pass = 0
    number_samples_warning = 0

    # Get MLST scheme to use
    scheme = 'unknown'
    species_genus, mlst_scheme_genus = None, None
    if not args.skipMLST and not args.skipSPAdes:
        scheme, species_genus, mlst_scheme_genus = mlst.getScheme(args.speciesExpected)
        # Print path to blastn
        mlst.getBlastPath()

    # Memory
    available_memory_GB = utils.get_free_memory() / (1024.0 ** 2)
    # Determine SPAdes maximum memory
    spadesMaxMemory = None
    if not args.skipSPAdes:
        print('')
        spadesMaxMemory = spades.define_memory(args.spadesMaxMemory, args.threads, available_memory_GB)
    # Determine .jar maximum memory
    jarMaxMemory = 'off'
    if not (args.skipTrimmomatic and (args.skipSPAdes or args.skipPilon)):
        print('')
        jarMaxMemory = utils.define_jar_max_memory(args.jarMaxMemory, args.threads, available_memory_GB)

    # Run INNUca for each sample
    sample_report_json = {}
    for sample in samples:
        sample_start_time = time.time()

        print('\n' + 'Sample: ' + sample + '\n')

        # Create sample outdir
        sample_outdir = os.path.abspath(os.path.join(outdir, sample, ''))
        if not os.path.isdir(sample_outdir):
            os.makedirs(sample_outdir)

        # Get fastq files
        fastq_files = utils.searchFastqFiles(os.path.join(inputDirectory, sample, ''), pairEnd_filesSeparation_list, False)
        if len(fastq_files) == 1:
            print('Only one fastq file was found: ' + str(fastq_files))
            print('Pair-End sequencing is required. Moving to the next sample')
            continue
        elif len(fastq_files) == 0:
            print('No compressed fastq files were found. Continue to the next sample')
            continue

        print('The following files will be used:')
        print(str(fastq_files) + '\n')

        # Run INNUca.py analysis
        run_successfully, pass_qc, run_report = \
            run_innuca(sample, sample_outdir, fastq_files, args, script_path, scheme, spadesMaxMemory,
                       jar_path_trimmomatic, jar_path_pilon, jarMaxMemory, trueCoverage_config, rematch_script,
                       species_genus, mlst_scheme_genus, spades_version=spades_version)

        # Save sample fail report
        utils.write_fail_report(os.path.join(sample_outdir, 'fail_report.txt'), run_report)

        # Save warning report
        write_warning_report(os.path.join(sample_outdir, 'warning_report.txt'), run_report)

        # Get raw reads files size
        fileSize = sum(os.path.getsize(fastq) for fastq in fastq_files)

        # Remove sample directory if it was created during the process
        if removeCreatedSamplesDirectories and not indir_same_outdir:
            utils.removeDirectory(os.path.join(inputDirectory, sample, ''))

        # Remove reads folder if --fastq was provided
        if args.fastq is not None:
            utils.removeDirectory(os.path.join(outdir, 'reads', ''))

        print('END ' + sample + ' analysis')
        time_taken = utils.runTime(sample_start_time)

        # Save run report
        warning, json_pass_qc = utils.write_sample_report(samples_report_path, sample, run_successfully, pass_qc, time_taken, fileSize, run_report)

        # Save runs statistics
        if run_successfully:
            number_samples_successfully += 1
        if pass_qc:
            if warning:
                number_samples_warning += 1
            else:
                number_samples_pass += 1

        sample_report_json[sample] = {'run_successfully': run_successfully, 'pass_qc': json_pass_qc, 'modules_run_report': run_report}

    # Save combine_samples_reports
    combine_reports.combine_reports(outdir, outdir, args.json, time_str, len(samples))

    # Save sample_report in json
    if args.json:
        import json
        with open(os.path.join(outdir, 'samples_report.' + time_str + '.json'), 'wt') as writer:
            json.dump(sample_report_json, writer)

    # Remove temporary folder with symlink to fastq files in case of --fastq use
    if args.inputDirectory is None and args.fastq is not None:
        utils.removeDirectory(os.path.join(inputDirectory, ''))

    # Run report
    print('\n' + 'END INNUca.py')
    print('\n' + 'Pipeline problems: {not_run_successfully}'
                 ' samples'.format(not_run_successfully=(len(samples) - number_samples_successfully)))
    print('\n' + 'FAIL: {number_samples_fail}'
                 ' samples'.format(number_samples_fail=(len(samples) - number_samples_pass - number_samples_warning)))
    print('\n' + 'WARNING: {number_samples_warning} samples'.format(number_samples_warning=number_samples_warning))
    print('\n' + 'PASS: {number_samples_pass} samples'.format(number_samples_pass=number_samples_pass))
    _ = utils.runTime(general_start_time)

    # Check whether INNUca.py run at least one sample successfully
    if number_samples_successfully == 0:
        sys.exit('No samples run successfully!')


def write_warning_report(warning_report_path, run_report):
    with open(warning_report_path, 'wt') as writer_warningReport:
        warnings = []
        for step in ('reads_Kraken', 'first_FastQC', 'Trimmomatic', 'second_FastQC', 'Pear', 'SPAdes',
                     'Assembly_Mapping', 'MLST', 'assembly_Kraken'):
            if len(run_report[step][4]) > 0:
                if step == 'first_FastQC' and \
                        run_report['second_FastQC'][1] is not False and \
                        len(run_report['second_FastQC'][4]) == 0:
                    continue
                else:
                    if len(run_report[step][4]) == 0:
                        continue
                    else:
                        warnings.append('#' + step)
                        for key, warning_reasons in run_report[step][4].items():
                            warnings.append('>' + str(key))
                            if isinstance(warning_reasons, (list, tuple)):
                                for reasons in warning_reasons:
                                    warnings.append(str(reasons))
                            else:
                                warnings.append(str(warning_reasons))
        writer_warningReport.write('\n'.join(warnings))


def get_sample_args_fastq(fastq_files_list, outdir, pairEnd_filesSeparation_list):
    new_indir = os.path.join(outdir, 'reads', '')
    utils.removeDirectory(new_indir)
    os.mkdir(new_indir)
    samples = []
    for fastq in fastq_files_list:
        fastq_link = os.path.join(new_indir, os.path.basename(fastq))
        os.symlink(fastq, fastq_link)
    samples, removeCreatedSamplesDirectories, indir_same_outdir = \
        utils.checkSetInputDirectory(new_indir, outdir, pairEnd_filesSeparation_list)
    return new_indir, samples, removeCreatedSamplesDirectories, indir_same_outdir


def get_samples(args_input_directory, args_fastq, outdir, pair_end_files_separation_list):
    samples, input_directory, remove_created_samples_directories, indir_same_outdir = None, None, None, None

    if args_fastq is None:
        # Check if input directory exists with fastq files and store samples name that have fastq files
        input_directory = os.path.abspath(os.path.join(args_input_directory, ''))
        print('')
        samples, remove_created_samples_directories, indir_same_outdir = \
            utils.checkSetInputDirectory(input_directory, outdir, pair_end_files_separation_list)
    elif args_input_directory is None:
        fastq_files = [os.path.abspath(fastq.name) for fastq in args_fastq]
        if fastq_files[0] == fastq_files[1]:
            sys.exit('Same fastq file provided twice')
        input_directory, samples, remove_created_samples_directories, indir_same_outdir = \
            get_sample_args_fastq(fastq_files, outdir, pair_end_files_separation_list)

    return samples, input_directory, remove_created_samples_directories, indir_same_outdir


def run_innuca(sample_name, outdir, fastq_files, args, script_path, scheme, spades_max_memory, jar_path_trimmomatic,
               jar_path_pilon, jar_max_memory, true_coverage_config, rematch_script, species_genus, mlst_scheme_genus,
               spades_version=None):
    threads = args.threads
    adapters_fasta = args.adapters
    if adapters_fasta is not None:
        adapters_fasta = os.path.abspath(adapters_fasta.name)
    genome_size = args.genomeSizeExpectedMb
    # run_successfully, pass_qc, time_taken, failing, warning, file_size
    skipped = [None, None, 0, {'sample': 'Skipped'}, {}, 'NA']
    not_run = [None, None, 0, {'sample': 'Not run'}, {}, 'NA']

    runs = {}

    # Run FastQ integrity check
    not_corruption_found, pass_qc, time_taken, failing, fastq_encoding, min_reads_length, max_reads_length = \
        fastQintegrity.runFastQintegrity(fastq_files, threads, outdir)
    runs['FastQ_Integrity'] = [not_corruption_found, pass_qc, time_taken, failing, {}, 'NA']

    pear_folder = None
    trimmomatic_folder = None
    if not_corruption_found:
        # Run Kraken
        # most_abundant_taxon_percent = None
        run_successfully_kraken = False
        run_successfully_estimated_coverage = False
        estimated_coverage = None
        run_successfully_true_coverage = False
        pass_qc_true_coverage = False

        trimmomatic_run_successfully = False
        if args.runKraken:
            print('\n'
                  '--runKraken set. Running Kraken for reads')
            run_successfully_kraken, pass_qc, time_taken, failing, warning, _ = \
                kraken(species=args.speciesExpected, files_to_classify=fastq_files, kraken_db=args.krakenDB,
                       files_type='fastq', outdir=outdir, version_kraken=version_kraken_global,
                       db_mem=args.krakenMemory, quick=args.krakenQuick, min_percent_covered=args.krakenMinCov,
                       max_unclassified_frag=args.krakenMaxUnclass, min_base_quality=args.krakenMinQual,
                       threads=threads)
            runs['reads_Kraken'] = [run_successfully_kraken, pass_qc, time_taken, failing, warning, 'NA']
        else:
            runs['reads_Kraken'] = skipped

        if args.runKraken and \
                (run_successfully_kraken and not pass_qc) and \
                not args.krakenProceed and \
                not args.krakenIgnoreQC:
            print('\n'
                  'This sample does not pass Kraken module QA/QC. It will not proceed with INNUca pipeline')
        else:
            # Run first Estimated Coverage
            if not args.skipEstimatedCoverage:
                # Check whether the Estimated Coverage output is already present
                report_file = os.path.join(outdir, 'coverage_report.txt')
                if os.path.isfile(report_file):
                    os.remove(report_file)
                # Run getEstimatedCoverage
                run_successfully_estimated_coverage, pass_qc, time_taken, failing, estimated_coverage = \
                    coverage.getEstimatedCoverage(fastq_files, genome_size, outdir, threads,
                                                  args.estimatedMinimumCoverage)
                runs['first_Coverage'] = [run_successfully_estimated_coverage, pass_qc, time_taken, failing, {}, 'NA']
            else:
                print('--skipEstimatedCoverage set. Skipping First Estimated Coverage analysis')
                runs['first_Coverage'] = skipped

            # # Correct first estimation coverage with Kraken percentage
            # # Does not seem to be a good idea (at least for Streptococcus agalactiae)
            # if args.runKraken and \
            #         (runs['Kraken'][0] and runs['Kraken'][1]) and \
            #         most_abundant_taxon_percent is not None and \
            #         estimated_coverage is not None:
            #     new_estimation = estimated_coverage * (most_abundant_taxon_percent / 100)
            #     print('\n'
            #           'Correct estimated coverage ({estimated}x) with Kraken taxon percentage'
            #           ' coverage ({percent}%): {new_estimation}x'.format(estimated=estimated_coverage,
            #                                                              percent=most_abundant_taxon_percent,
            #                                                              new_estimation=new_estimation))
            #     estimated_coverage = new_estimation

            if args.skipEstimatedCoverage or (run_successfully_estimated_coverage and
                                              not estimated_coverage < args.estimatedMinimumCoverage):
                if not args.skipTrueCoverage and true_coverage_config is not None:
                    # Run True Coverage
                    run_successfully_true_coverage, pass_qc_true_coverage, time_taken, failing, _ = \
                        trueCoverage.run_true_coverage(sample_name, fastq_files, true_coverage_config['reference_file'],
                                                       threads, outdir,
                                                       true_coverage_config['length_extra_seq'],
                                                       true_coverage_config['minimum_depth_presence'],
                                                       true_coverage_config['minimum_depth_call'],
                                                       true_coverage_config['minimum_depth_frequency_dominant_allele'],
                                                       true_coverage_config['minimum_gene_coverage'], False,
                                                       true_coverage_config['minimum_gene_identity'],
                                                       true_coverage_config, rematch_script, num_map_loc=1,
                                                       bowtie_algorithm=args.trueCoverageBowtieAlgo,
                                                       clean_run_rematch=True)
                    runs['trueCoverage_ReMatCh'] = [run_successfully_true_coverage, pass_qc_true_coverage, time_taken,
                                                    failing, {}, 'NA']
                else:
                    print('\n' + '--skipTrueCoverage set. Skipping True coverage analysis')
                    runs['trueCoverage_ReMatCh'] = skipped

                if args.skipTrueCoverage or true_coverage_config is None or args.trueCoverageProceed or \
                        (run_successfully_true_coverage and pass_qc_true_coverage):
                    # Run first FastQC
                    nts2clip_based_nts_content = None
                    if not args.skipFastQC:
                        run_successfully, pass_qc, time_taken, failing, warning, maximum_reads_length, \
                            nts2clip_based_nts_content = fastqc.runFastQCanalysis(outdir, threads, adapters_fasta,
                                                                                  fastq_files, args.fastQCkeepFiles,
                                                                                  'first_run')
                        runs['first_FastQC'] = [run_successfully, pass_qc, time_taken, failing, warning, 'NA']
                    else:
                        print('--skipFastQC set. Skipping First FastQC analysis')
                        runs['first_FastQC'] = skipped

                    # Run Trimmomatic
                    not_empty_fastq = True
                    if not args.skipTrimmomatic:
                        run_successfully, not_empty_fastq, time_taken, failing, paired_reads, trimmomatic_folder, \
                            file_size, warning = trimmomatic.runTrimmomatic(jar_path_trimmomatic, sample_name, outdir,
                                                                            threads, adapters_fasta, script_path,
                                                                            args.doNotSearchAdapters, fastq_files,
                                                                            max_reads_length, args.doNotTrimCrops,
                                                                            args.trimCrop, args.trimHeadCrop,
                                                                            args.trimLeading, args.trimTrailing,
                                                                            args.trimSlidingWindow, args.trimMinLength,
                                                                            nts2clip_based_nts_content, jar_max_memory,
                                                                            fastq_encoding)
                        runs['Trimmomatic'] = [run_successfully, None, time_taken, failing, warning, file_size]
                        trimmomatic_run_successfully = run_successfully

                        if run_successfully and not_empty_fastq:
                            fastq_files = paired_reads
                            min_reads_length = args.trimMinLength

                            # Run second Estimated Coverage
                            if not args.skipEstimatedCoverage:
                                run_successfully_estimated_coverage, pass_qc, time_run, failing, estimated_coverage = \
                                    coverage.getEstimatedCoverage(fastq_files, genome_size, outdir, threads,
                                                                  args.estimatedMinimumCoverage)
                                runs['second_Coverage'] = [run_successfully_estimated_coverage, pass_qc, time_run,
                                                           failing, {}, 'NA']
                            else:
                                print('--skipEstimatedCoverage set. Skipping Second Estimated Coverage analysis')
                                runs['second_Coverage'] = skipped

                            if args.skipEstimatedCoverage or (run_successfully_estimated_coverage and
                                                              not estimated_coverage < args.estimatedMinimumCoverage):
                                # Run second FastQC
                                if not args.skipFastQC:
                                    run_successfully, pass_qc, time_taken, failing, warning, maximum_reads_length, \
                                        nts2clip_based_nts_content = fastqc.runFastQCanalysis(outdir, threads,
                                                                                              adapters_fasta,
                                                                                              fastq_files,
                                                                                              args.fastQCkeepFiles,
                                                                                              'second_run')
                                    runs['second_FastQC'] = [run_successfully, pass_qc, time_taken, failing, warning,
                                                             'NA']
                                    if run_successfully:
                                        max_reads_length = maximum_reads_length
                                else:
                                    print('--skipFastQC set. Skipping Second FastQC analysis')
                                    runs['second_FastQC'] = skipped
                            else:
                                print('\n'
                                      'Estimated coverage is too lower (< {estimatedMinimumCoverage}x). This sample'
                                      ' will not proceed with INNUca'
                                      ' pipeline'.format(estimatedMinimumCoverage=args.estimatedMinimumCoverage))
                                runs['second_FastQC'] = skipped
                        else:
                            print('Trimmomatic did not run successfully or return zero reads! Skipping Second Estimated'
                                  ' Coverage analysis and FastQC analysis')
                            runs['second_Coverage'] = skipped
                            runs['second_FastQC'] = skipped
                    else:
                        print('--skipTrimmomatic set. Skipping Trimmomatic, but also Second FastQC analysis and Second'
                              ' Estimated Coverage analysis')
                        runs['Trimmomatic'] = skipped
                        runs['second_Coverage'] = skipped
                        runs['second_FastQC'] = skipped

                    if not args.skipFastQC and \
                            (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and
                                                          runs['first_FastQC'][1])) is False and \
                            not not_empty_fastq and not args.fastQCproceed:
                        print('\n'
                              'This sample does not pass FastQC module QA/QC. It will not proceed with INNUca pipeline')
                else:
                    print('\n'
                          'This sample does not pass True Coverage module QA/QC. This sample will not proceed with'
                          ' INNUca pipeline')
            else:
                print('\n'
                      'Estimated coverage is too lower (< {estimatedMinimumCoverage}x). This sample will not proceed'
                      ' with INNUca pipeline'.format(estimatedMinimumCoverage=args.estimatedMinimumCoverage))

        continue_second_part = False
        if not args.runKraken or \
                (runs['reads_Kraken'][0] is True and runs['reads_Kraken'][1] is True) or \
                args.krakenProceed or \
                args.krakenIgnoreQC:
            if args.skipEstimatedCoverage or (run_successfully_estimated_coverage and
                                              not estimated_coverage < args.estimatedMinimumCoverage):
                if args.skipTrueCoverage or true_coverage_config is None or args.trueCoverageProceed or \
                        (run_successfully_true_coverage and pass_qc_true_coverage):
                    if args.skipFastQC or (runs['second_FastQC'][1] or
                                           (runs['second_FastQC'][1] is None and
                                            runs['first_FastQC'][1])) is not False or \
                            args.fastQCproceed:
                        continue_second_part = True

        if continue_second_part:
            unassembled_pe_reads = None
            assembled_se_reads = None
            # Run Pear
            if args.runPear:
                print('--runPear set. Running Pear')
                pear_min_overlap = pear.determine_minimum_overlap(args.pearMinOverlap, min_reads_length,
                                                                  max_reads_length)
                run_successfully, pass_qc, time_taken, failing, unassembled_pe_reads, assembled_se_reads, \
                    pear_folder, warning = pear.runPear(fastq_files, threads, outdir, sample_name,
                                                        fastq_encoding, trimmomatic_run_successfully,
                                                        pear_min_overlap)
                runs['Pear'] = [run_successfully, pass_qc, time_taken, failing, warning, 'NA']
            else:
                runs['Pear'] = skipped

            # Run SPAdes
            if not args.skipSPAdes:
                run_successfully, pass_qc, time_taken, failing, contigs_spades, warning = \
                    spades.run_spades(sample_name, outdir, threads,
                                      unassembled_pe_reads if unassembled_pe_reads is not None else fastq_files,
                                      args.spadesNotUseCareful, spades_max_memory,
                                      args.spadesMinCoverageAssembly, args.spadesMinContigsLength, genome_size,
                                      args.spadesKmers, max_reads_length, args.spadesDefaultKmers,
                                      args.spadesMinKmerCovContigs, assembled_se_reads, args.saveExcludedContigs,
                                      args.maxNumberContigs, args.keepSPAdesScaffolds, spades_version=spades_version,
                                      estimated_coverage=estimated_coverage,
                                      spades_not_use_isolate=args.spadesNotUseIsolate)
                runs['SPAdes'] = [run_successfully, pass_qc, time_taken, failing, warning, 'NA']

                if run_successfully:
                    contigs = contigs_spades

                    # Run Assembly Mapping check
                    bam_file = None
                    original_bam = None
                    assembly_mapping_folder = None
                    possible_assemblies_bam_remove = {}
                    if not args.skipAssemblyMapping:
                        run_successfully, pass_qc, time_taken, failing, assembly_filtered, bam_file, \
                            assembly_mapping_folder, warning, original_bam = \
                            assembly_mapping.run_assembly_mapping(fastq_files=fastq_files, reference_file=contigs,
                                                                  outdir=outdir, estimated_genome_size_mb=genome_size,
                                                                  max_number_contigs=args.maxNumberContigs,
                                                                  save_excluded_contigs=args.saveExcludedContigs,
                                                                  min_coverage_assembly=args.assemblyMinCoverageContigs,
                                                                  keep_bam=args.keepBAM, threads=threads)
                        runs['Assembly_Mapping'] = [run_successfully, pass_qc, time_taken, failing, warning,
                                                    'NA']

                        if run_successfully:
                            # Assembly to remove
                            if not args.keepIntermediateAssemblies:
                                if os.path.isfile(contigs_spades) and \
                                        assembly_filtered is not None and \
                                        assembly_filtered != contigs_spades:
                                    if not args.keepBAM:
                                        os.remove(contigs_spades)
                                    else:
                                        possible_assemblies_bam_remove['assembly_mapping'] = contigs_spades

                            if assembly_filtered is not None and \
                                    assembly_filtered != contigs_spades and \
                                    os.path.isfile(assembly_filtered):
                                contigs = assembly_filtered
                    else:
                        print('--skipAssemblyMapping set. Skipping Assembly Mapping check')
                        runs['Assembly_Mapping'] = skipped

                    # Run Pilon
                    pilon_new_bam = False
                    pilon_bam = None
                    if not args.skipPilon:
                        run_successfully, _, time_taken, failing, assembly_polished, pilon_folder, pilon_new_bam, \
                            pilon_bam = pilon.run_pilon(jar_path_pilon=jar_path_pilon, assembly=contigs,
                                                        fastq_files=fastq_files, outdir=outdir,
                                                        jar_max_memory=jar_max_memory, alignment_file=bam_file,
                                                        keep_bam=args.keepBAM, threads=threads)
                        runs['Pilon'] = [run_successfully, None, time_taken, failing, {}, 'NA']

                        if run_successfully:
                            if not args.keepIntermediateAssemblies:
                                if os.path.isfile(contigs) and \
                                        assembly_polished is not None and \
                                        os.path.isfile(assembly_polished):
                                    if not args.keepBAM:
                                        os.remove(contigs)
                                    else:
                                        if not pilon_new_bam:
                                            possible_assemblies_bam_remove['pilon'] = contigs

                            if assembly_polished is not None and \
                                    os.path.isfile(assembly_polished):
                                contigs = assembly_polished

                        if not args.pilonKeepFiles and os.path.isdir(pilon_folder):
                            utils.removeDirectory(pilon_folder)

                    else:
                        print('--skipPilon set. Skipping Pilon correction')
                        runs['Pilon'] = skipped

                    if not args.keepBAM:
                        if bam_file is not None:
                            if os.path.isfile(bam_file):
                                os.remove(bam_file)
                            if os.path.isfile(bam_file + '.bai'):
                                os.remove(bam_file + '.bai')

                        if original_bam is not None and os.path.isfile(original_bam):
                            os.remove(original_bam)

                        if pilon_bam is not None and os.path.isfile(pilon_bam):
                            os.remove(pilon_bam)

                        if 'assembly_mapping' in possible_assemblies_bam_remove and \
                                os.path.isfile(possible_assemblies_bam_remove['assembly_mapping']):
                            os.remove(possible_assemblies_bam_remove['assembly_mapping'])
                        if 'pilon' in possible_assemblies_bam_remove and \
                                os.path.isfile(possible_assemblies_bam_remove['pilon']):
                            os.remove(possible_assemblies_bam_remove['pilon'])
                    else:
                        if pilon_new_bam:
                            if bam_file is not None:
                                if os.path.isfile(bam_file):
                                    os.remove(bam_file)
                                if os.path.isfile(bam_file + '.bai'):
                                    os.remove(bam_file + '.bai')

                            if original_bam is not None and os.path.isfile(original_bam):
                                os.remove(original_bam)

                            if 'assembly_mapping' in possible_assemblies_bam_remove and \
                                    os.path.isfile(possible_assemblies_bam_remove['assembly_mapping']):
                                os.remove(possible_assemblies_bam_remove['assembly_mapping'])
                        else:
                            if original_bam is not None and os.path.isfile(original_bam) and \
                                    bam_file is not None and os.path.isfile(bam_file):
                                os.remove(bam_file)
                            if 'pilon' in possible_assemblies_bam_remove and \
                                    os.path.isfile(possible_assemblies_bam_remove['pilon']):
                                os.remove(possible_assemblies_bam_remove['pilon'])

                    if not args.skipAssemblyMapping:
                        utils.removeDirectory(assembly_mapping_folder)

                    print('\n' + 'Final assembly: ' + contigs)
                    with open(os.path.join(outdir, 'final_assembly.txt'), 'wt') as writer:
                        writer.write(contigs + '\n')

                    # Run MLST
                    if not args.skipMLST:
                        run_successfully, pass_qc, time_taken, failing, warning = \
                            mlst.runMlst(contigs, scheme, outdir, species_genus, mlst_scheme_genus)
                        runs['MLST'] = [run_successfully, pass_qc, time_taken, failing, warning, 'NA']
                    else:
                        print('--skipMLST set. Skipping MLST analysis')
                        runs['MLST'] = skipped

                    # Run Kraken
                    if args.runKraken:
                        print('\n'
                              '--runKraken set. Running Kraken for assembly')
                        run_successfully, pass_qc, time_taken, failing, warning, _ = \
                            kraken(species=args.speciesExpected, files_to_classify=[contigs], kraken_db=args.krakenDB,
                                   files_type='fasta', outdir=outdir, version_kraken=version_kraken_global,
                                   db_mem=args.krakenMemory, quick=args.krakenQuick,
                                   min_percent_covered=args.krakenMinCov,
                                   max_unclassified_frag=args.krakenMaxUnclass, min_base_quality=args.krakenMinQual,
                                   threads=threads)
                        runs['assembly_Kraken'] = [run_successfully, pass_qc, time_taken, failing, warning, 'NA']
                    else:
                        runs['assembly_Kraken'] = skipped

                    # Run insert_size
                    if args.runInsertSize:
                        print('\n'
                              '--runInsertSize set. Running insert_size')
                        run_successfully, _, time_taken, failing = \
                            insert_size(sample_name=sample_name, reference=contigs,
                                        fastq=fastq_files, outdir=outdir, threads=threads, dist=args.insertSizeDist)
                        runs['insert_size'] = [run_successfully, None, time_taken, failing, {}, 'NA']
                    else:
                        runs['insert_size'] = skipped
                else:
                    print('SPAdes did not run successfully! Skipping Pilon correction, Assembly Mapping check,'
                          ' MLST and Kraken (assembly) analysis and insert size determination')
            else:
                print('--skipSPAdes set. Skipping SPAdes, Pilon correction, Assembly Mapping check and MLST and Kraken'
                      ' (assembly) analysis and insert size determination')
                runs['SPAdes'] = skipped
                runs['Assembly_Mapping'] = skipped
                runs['Pilon'] = skipped
                runs['MLST'] = skipped
                runs['assembly_Kraken'] = skipped
                runs['insert_size'] = skipped
    else:
        print('Moving to the next sample')

    for step in ('reads_Kraken', 'first_Coverage', 'trueCoverage_ReMatCh', 'first_FastQC', 'Trimmomatic',
                 'second_Coverage', 'second_FastQC', 'Pear', 'SPAdes', 'Assembly_Mapping', 'Pilon', 'MLST',
                 'assembly_Kraken', 'insert_size'):
        if step not in runs:
            runs[step] = not_run

    # Remove Pear directory
    if not args.pearKeepFiles and pear_folder is not None:
        utils.removeDirectory(pear_folder)
    # Remove Trimmomatic directory with cleaned reads
    if not args.trimKeepFiles and trimmomatic_folder is not None:
        utils.removeDirectory(trimmomatic_folder)

    # Check run
    run_successfully = all(runs[step][0] or runs[step][0] is None for step in runs)

    pass_fastq_integrity = runs['FastQ_Integrity'][0]
    pass_reads_kraken = runs['reads_Kraken'][1] is not False or args.krakenIgnoreQC
    pass_cov = (runs['second_Coverage'][1] or (runs['second_Coverage'][1] is None and
                                               runs['first_Coverage'][1])) is not False
    pass_true_cov = runs['trueCoverage_ReMatCh'][1] is not False or args.trueCoverageIgnoreQC
    pass_fastqc = (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and
                                                runs['first_FastQC'][1])) is not False
    # pass_trimmomatic = runs['Trimmomatic'][1] is not False
    # pass_pear = runs['Pear'][1] is not False
    # pass_spades = runs['SPAdes'][1] is not False or runs['Assembly_Mapping'][1] is True
    pass_spades = runs['SPAdes'][1] is not False
    pass_assembly_mapping = runs['Assembly_Mapping'][1] is not False
    pass_pilon = runs['Pilon'][0] is not False
    pass_mlst = runs['MLST'][1] is not False or args.mlstIgnoreQC
    pass_assembly_kraken = runs['assembly_Kraken'][1] is not False or args.krakenIgnoreQC
    pass_qc = all([pass_fastq_integrity, pass_reads_kraken, pass_cov, pass_true_cov, pass_fastqc, pass_spades,
                   pass_assembly_mapping, pass_pilon, pass_mlst, pass_assembly_kraken])

    return run_successfully, pass_qc, runs


if __name__ == "__main__":
    main()
