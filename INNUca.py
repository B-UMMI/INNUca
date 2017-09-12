#!/usr/bin/env python

# -*- coding: utf-8 -*-


"""
INNUca - Reads Control and Assembly
INNUca.py - INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search
<https://github.com/B-UMMI/INNUca>

Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: June 21, 2017

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
import modules.fastQintegrity as fastQintegrity
import modules.estimated_coverage as coverage
import modules.fastqc as fastqc
import modules.trimmomatic as trimmomatic
import modules.pear as pear
import modules.spades as spades
import modules.assembly_mapping as assembly_mapping
import modules.pilon as pilon
import modules.mlst as mlst
import modules.trueCoverage_rematch as trueCoverage
import modules.combine_reports as combine_reports
import time
import os
import sys


def get_trueCoverage_config(skipTrueCoverage, trueConfigFile, speciesExpected, script_path):
    trueCoverage_config = None
    if not skipTrueCoverage:
        trueCoverage_reference = None
        trueCoverage_config_file = None
        trueCoverage_config = None

        if trueConfigFile is None:
            print 'No trueCoverage_ReMatCh config file was provided. Search for default files'
            trueCoverage_config_file, trueCoverage_reference = trueCoverage.check_existing_default_config(speciesExpected, script_path)
        else:
            trueCoverage_config_file = trueConfigFile

        if trueCoverage_config_file is not None:
            trueCoverage_config = trueCoverage.parse_config(trueCoverage_config_file)
        if trueConfigFile is None and trueCoverage_config is not None:
            trueCoverage_config['reference_file'] = trueCoverage_reference

        if trueCoverage_config is not None:
            print 'The following trueCoverage_ReMatCh config file will be used: ' + trueCoverage_config_file
            print 'The following trueCoverage_ReMatCh reference file will be used: ' + trueCoverage_config['reference_file'] + '\n'
        else:
            print 'No trueCoverage_ReMatCh config file was found'
    return trueCoverage_config


def include_rematch_dependencies_path(doNotUseProvidedSoftware):
    command = ['which', 'rematch.py']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    if run_successfully:
        rematch = stdout.splitlines()[0]
        path_variable = os.environ['PATH']
        script_folder = os.path.dirname(rematch)
        if not doNotUseProvidedSoftware:
            bcftools = os.path.join(script_folder, 'src', 'bcftools-1.3.1', 'bin')
            os.environ['PATH'] = str(':'.join([bcftools, path_variable]))


def main():
    version = '3.1'
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

    print '\n' + '==========> INNUca.py <=========='
    print '\n' + 'Program start: ' + time.ctime()

    # Tells where the logfile will be stored
    if not args.noLog:
        print '\n' + 'LOGFILE:'
        print sys.stdout.getLogFile()

    # Print command
    print '\n' + 'COMMAND:'
    script_path = os.path.abspath(sys.argv[0])
    print sys.executable + ' ' + script_path + ' ' + ' '.join(sys.argv[1:])

    # Print directory where programme was lunch
    print '\n' + 'PRESENT DIRECTORY:'
    print os.getcwd()

    # Print program version
    print '\n' + 'VERSION INNUca.py:'
    utils.scriptVersionGit(version, os.getcwd(), script_path, args.noGitInfo)

    # Get CPU information
    utils.get_cpu_information(outdir, time_str)

    # Get trueCoverage_ReMatCh settings
    trueCoverage_config = get_trueCoverage_config(args.skipTrueCoverage, args.trueConfigFile.name if args.trueConfigFile is not None else None, args.speciesExpected, script_path)

    # Check programms
    programs_version_dictionary = {}
    programs_version_dictionary['gunzip'] = ['--version', '>=', '1.6']

    # Java check first for java dependents check next
    if not (args.skipFastQC and args.skipTrimmomatic and (args.skipPilon or args.skipSPAdes)):
        # programs_version_dictionary['java'] = ['-version', '>=', '1.8']
        programs_version_dictionary['java'] = [None, '>=', '1.8']  # For OpenJDK compatibility
    missingPrograms, programs_version_dictionary = utils.checkPrograms(programs_version_dictionary)
    if len(missingPrograms) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))

    if not args.skipTrueCoverage or trueCoverage_config is not None:
        include_rematch_dependencies_path(args.doNotUseProvidedSoftware)
        programs_version_dictionary['rematch.py'] = ['--version', '>=', '3.2']
        programs_version_dictionary['bcftools'] = ['--version', '==', '1.3.1']
    if not (args.skipTrueCoverage and ((args.skipAssemblyMapping and args.skipPilon) or args.skipSPAdes)):
        programs_version_dictionary['bowtie2'] = ['--version', '>=', '2.2.9']
        programs_version_dictionary['samtools'] = ['--version', '==', '1.3.1']
    if not args.skipFastQC:
        programs_version_dictionary['fastqc'] = ['--version', '==', '0.11.5']
    if not args.skipTrimmomatic:
        programs_version_dictionary['trimmomatic-0.36.jar'] = ['-version', '==', '0.36']
    if args.runPear:
        programs_version_dictionary['pear'] = ['--version', '>=', '0.9.10']
    if not args.skipSPAdes:
        programs_version_dictionary['spades.py'] = ['--version', '>=', '3.9.0']
    if not (args.skipPilon or args.skipSPAdes):
        programs_version_dictionary['pilon-1.18.jar'] = ['--version', '==', '1.18']
    if not (args.skipMLST or args.skipSPAdes):
        programs_version_dictionary['mlst'] = ['--version', '>=', '2.4']

    # Set and print PATH variable
    utils.setPATHvariable(args, script_path)

    missingPrograms, programs_version_dictionary = utils.checkPrograms(programs_version_dictionary)
    if len(missingPrograms) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))

    # .jar paths
    jar_path_trimmomatic = None
    if not args.skipTrimmomatic:
        jar_path_trimmomatic = programs_version_dictionary['trimmomatic-0.36.jar'][3]

    jar_path_pilon = None
    if not args.skipPilon and not args.skipSPAdes:
        jar_path_pilon = programs_version_dictionary['pilon-1.18.jar'][3]

    rematch_script = None
    # ReMatCh path
    if not args.skipTrueCoverage:
        rematch_script = programs_version_dictionary['rematch.py'][3]

    # pairEnd_filesSeparation_list = args.pairEnd_filesSeparation
    pairEnd_filesSeparation_list = None
    samples, inputDirectory, removeCreatedSamplesDirectories, indir_same_outdir = get_samples(args.inputDirectory, args.fastq, outdir, pairEnd_filesSeparation_list)

    # Start running the analysis
    print '\n' + 'RUNNING INNUca.py'

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
        print ''
        spadesMaxMemory = spades.define_memory(args.spadesMaxMemory, args.threads, available_memory_GB)
    # Determine .jar maximum memory
    jarMaxMemory = 'off'
    if not (args.skipTrimmomatic and (args.skipSPAdes or args.skipPilon)):
        print ''
        jarMaxMemory = utils.define_jar_max_memory(args.jarMaxMemory, args.threads, available_memory_GB)

    # Run INNUca for each sample
    sample_report_json = {}
    for sample in samples:
        sample_start_time = time.time()

        print '\n' + 'Sample: ' + sample + '\n'

        # Create sample outdir
        sample_outdir = os.path.abspath(os.path.join(outdir, sample, ''))
        if not os.path.isdir(sample_outdir):
            os.makedirs(sample_outdir)

        # Get fastq files
        fastq_files = utils.searchFastqFiles(os.path.join(inputDirectory, sample, ''), pairEnd_filesSeparation_list, False)
        if len(fastq_files) == 1:
            print 'Only one fastq file was found: ' + str(fastq_files)
            print 'Pair-End sequencing is required. Moving to the next sample'
            continue
        elif len(fastq_files) == 0:
            print 'No compressed fastq files were found. Continue to the next sample'
            continue

        print 'The following files will be used:'
        print str(fastq_files) + '\n'

        # Run INNUca.py analysis
        run_successfully, pass_qc, run_report = run_INNUca(sample, sample_outdir, fastq_files, args, script_path, scheme, spadesMaxMemory, jar_path_trimmomatic, jar_path_pilon, jarMaxMemory, trueCoverage_config, rematch_script, species_genus, mlst_scheme_genus)

        # Save sample fail report
        utils.write_fail_report(os.path.join(sample_outdir, 'fail_report.txt'), run_report)

        # Save warning report
        write_warning_report(os.path.join(sample_outdir, 'warning_report.txt'), run_report)

        # Get raw reads files size
        fileSize = sum(os.path.getsize(fastq) for fastq in fastq_files)

        # Remove sample directory if it was created during the process
        if removeCreatedSamplesDirectories and not indir_same_outdir:
            utils.removeDirectory(os.path.join(inputDirectory, sample, ''))

        print 'END ' + sample + ' analysis'
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
    print '\n' + 'END INNUca.py'
    print '\n' + 'Pipeline problems: {not_run_successfully} samples'.format(not_run_successfully=(len(samples) - number_samples_successfully))
    print '\n' + 'FAIL: {number_samples_fail} samples'.format(number_samples_fail=(len(samples) - number_samples_pass - number_samples_warning))
    print '\n' + 'WARNING: {number_samples_warning} samples'.format(number_samples_warning=number_samples_warning)
    print '\n' + 'PASS: {number_samples_pass} samples'.format(number_samples_pass=number_samples_pass)
    time_taken = utils.runTime(general_start_time)
    del time_taken

    # Check whether INNUca.py run at least one sample successfully
    if number_samples_successfully == 0:
        sys.exit('No samples run successfully!')


def write_warning_report(warning_report_path, run_report):
    with open(warning_report_path, 'wt') as writer_warningReport:
        warnings = []
        for step in ('first_FastQC', 'second_FastQC', 'Pear', 'SPAdes', 'Assembly_Mapping', 'MLST'):
            if len(run_report[step][4]) > 0:
                if step == 'first_FastQC' and run_report['second_FastQC'][1] is not False and len(run_report['second_FastQC'][4]) == 0:
                    continue
                else:
                    if run_report[step][4] == 'NA':
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
    samples, removeCreatedSamplesDirectories, indir_same_outdir = utils.checkSetInputDirectory(new_indir, outdir, pairEnd_filesSeparation_list)
    return new_indir, samples, removeCreatedSamplesDirectories, indir_same_outdir


def get_samples(args_inputDirectory, args_fastq, outdir, pairEnd_filesSeparation_list):
    if args_fastq is None:
        # Check if input directory exists with fastq files and store samples name that have fastq files
        inputDirectory = os.path.abspath(os.path.join(args_inputDirectory, ''))
        print ''
        samples, removeCreatedSamplesDirectories, indir_same_outdir = utils.checkSetInputDirectory(inputDirectory, outdir, pairEnd_filesSeparation_list)
    elif args_inputDirectory is None:
        fastq_files = [os.path.abspath(fastq.name) for fastq in args_fastq]
        if fastq_files[0] == fastq_files[1]:
            sys.exit('Same fastq file provided twice')
        inputDirectory, samples, removeCreatedSamplesDirectories, indir_same_outdir = get_sample_args_fastq(fastq_files, outdir, pairEnd_filesSeparation_list)

    return samples, inputDirectory, removeCreatedSamplesDirectories, indir_same_outdir


def run_INNUca(sampleName, outdir, fastq_files, args, script_path, scheme, spadesMaxMemory, jar_path_trimmomatic, jar_path_pilon, jarMaxMemory, trueCoverage_config, rematch_script, species_genus, mlst_scheme_genus):
    threads = args.threads
    adaptersFasta = args.adapters
    if adaptersFasta is not None:
        adaptersFasta = os.path.abspath(adaptersFasta.name)
    genomeSize = args.genomeSizeExpectedMb
    skipped = [None, None, 0, {'sample': 'Skipped'}]
    not_run = [None, None, 0, {'sample': 'Not run'}]

    runs = {}

    # Run FastQ integrity check
    not_corruption_found, pass_qc, time_taken, failing, fastq_encoding, min_reads_length, max_reads_length = fastQintegrity.runFastQintegrity(fastq_files, threads, outdir)
    runs['FastQ_Integrity'] = [not_corruption_found, pass_qc, time_taken, failing]

    if not_corruption_found:
        # Run first Estimated Coverage
        run_successfully_estimatedCoverage = False
        estimatedCoverage = None
        run_successfully_trueCoverage = False
        pass_qc_trueCoverage = False
        if not args.skipEstimatedCoverage:
            # Check whether the Estimated Coverage output is already present
            report_file = os.path.join(outdir, 'coverage_report.txt')
            if os.path.isfile(report_file):
                os.remove(report_file)
            # Run getEstimatedCoverage
            run_successfully_estimatedCoverage, pass_qc, time_taken, failing, estimatedCoverage = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir, threads, args.estimatedMinimumCoverage)
            runs['first_Coverage'] = [run_successfully_estimatedCoverage, pass_qc, time_taken, failing]
        else:
            print '--skipEstimatedCoverage set. Skipping First Estimated Coverage analysis'
            runs['first_Coverage'] = skipped

        trimmomatic_run_successfully = False

        if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
            if not args.skipTrueCoverage and trueCoverage_config is not None:
                # Run True Coverage
                run_successfully_trueCoverage, pass_qc_trueCoverage, time_taken, failing = trueCoverage.runTrueCoverage(sampleName, fastq_files, trueCoverage_config['reference_file'], threads, outdir, trueCoverage_config['length_extra_seq'], trueCoverage_config['minimum_depth_presence'], trueCoverage_config['minimum_depth_call'], trueCoverage_config['minimum_depth_frequency_dominant_allele'], trueCoverage_config['minimum_gene_coverage'], False, False, 1, trueCoverage_config['minimum_gene_identity'], trueCoverage_config, rematch_script)
                runs['trueCoverage_ReMatCh'] = [run_successfully_trueCoverage, pass_qc_trueCoverage, time_taken, failing]
            else:
                print '\n' + '--skipTrueCoverage set. Skipping True coverage analysis'
                runs['trueCoverage_ReMatCh'] = skipped

            if args.skipTrueCoverage or trueCoverage_config is None or (run_successfully_trueCoverage and pass_qc_trueCoverage):
                # Run first FastQC
                nts2clip_based_ntsContent = None
                if not args.skipFastQC:
                    run_successfully, pass_qc, time_taken, failing, warning, maximum_reads_length, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files, args.fastQCkeepFiles, 'first_run')
                    runs['first_FastQC'] = [run_successfully, pass_qc, time_taken, failing, warning]
                else:
                    print '--skipFastQC set. Skipping First FastQC analysis'
                    runs['first_FastQC'] = skipped + ['NA']

                # Run Trimmomatic
                if not args.skipTrimmomatic:
                    run_successfully, not_empty_fastq, time_taken, failing, paired_reads, trimmomatic_folder, fileSize = trimmomatic.runTrimmomatic(jar_path_trimmomatic, sampleName, outdir, threads, adaptersFasta, script_path, args.doNotSearchAdapters, fastq_files, max_reads_length, args.doNotTrimCrops, args.trimCrop, args.trimHeadCrop, args.trimLeading, args.trimTrailing, args.trimSlidingWindow, args.trimMinLength, nts2clip_based_ntsContent, jarMaxMemory, fastq_encoding)
                    runs['Trimmomatic'] = [run_successfully, None, time_taken, failing, fileSize]
                    trimmomatic_run_successfully = run_successfully

                    if run_successfully and not_empty_fastq:
                        fastq_files = paired_reads
                        min_reads_length = args.trimMinLength

                        # Run second Estimated Coverage
                        if not args.skipEstimatedCoverage:
                            run_successfully_estimatedCoverage, pass_qc, time_taken, failing, estimatedCoverage = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir, threads, args.estimatedMinimumCoverage)
                            runs['second_Coverage'] = [run_successfully_estimatedCoverage, pass_qc, time_taken, failing]
                        else:
                            print '--skipEstimatedCoverage set. Skipping Second Estimated Coverage analysis'
                            runs['second_Coverage'] = skipped

                        if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
                            # Run second FastQC
                            if not args.skipFastQC:
                                run_successfully, pass_qc, time_taken, failing, warning, maximum_reads_length, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files, args.fastQCkeepFiles, 'second_run')
                                runs['second_FastQC'] = [run_successfully, pass_qc, time_taken, failing, warning]
                                if run_successfully:
                                    max_reads_length = maximum_reads_length
                            else:
                                print '--skipFastQC set. Skipping Second FastQC analysis'
                                runs['second_FastQC'] = skipped + ['NA']
                        else:
                            print '\n' + 'Estimated coverage is too lower (< ' + str(args.estimatedMinimumCoverage) + 'x). This sample will not proceed with INNUca pipeline'
                            runs['second_FastQC'] = not_run + ['NA']
                            runs['Pear'] = not_run + ['NA']
                            runs['SPAdes'] = not_run + ['NA']
                            runs['Assembly_Mapping'] = not_run + ['NA']
                            runs['Pilon'] = not_run
                            runs['MLST'] = not_run + ['NA']
                    else:
                        print 'Trimmomatic did not run successfully or return zero reads! Skipping Second Estimated Coverage analysis and FastQC analysis'
                        runs['second_Coverage'] = skipped
                        runs['second_FastQC'] = skipped + ['NA']

                else:
                    print '--skipTrimmomatic set. Skipping Trimmomatic, but also Second FastQC analysis and Second Estimated Coverage analysis'
                    runs['Trimmomatic'] = skipped + ['NA']
                    runs['second_Coverage'] = skipped
                    runs['second_FastQC'] = skipped + ['NA']

                if not args.skipFastQC and (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and runs['first_FastQC'][1])) is False and not args.fastQCproceed:
                    print '\n' + 'This sample does not pass FastQC module QA/QC. It will not proceed with INNUca pipeline'
                    runs['Pear'] = not_run + ['NA']
                    runs['SPAdes'] = not_run + ['NA']
                    runs['Assembly_Mapping'] = not_run + ['NA']
                    runs['Pilon'] = not_run
                    runs['MLST'] = not_run + ['NA']
            else:
                print '\n' + 'This sample does not pass True Coverage module QA/QC. This sample will not proceed with INNUca pipeline'
                runs['first_FastQC'] = not_run + ['NA']
                runs['Trimmomatic'] = not_run + ['NA']
                runs['second_Coverage'] = not_run
                runs['second_FastQC'] = not_run + ['NA']
                runs['Pear'] = not_run + ['NA']
                runs['SPAdes'] = not_run + ['NA']
                runs['Assembly_Mapping'] = not_run + ['NA']
                runs['Pilon'] = not_run
                runs['MLST'] = not_run + ['NA']

        else:
            print '\n' + 'Estimated coverage is too lower (< ' + str(args.estimatedMinimumCoverage) + 'x). This sample will not proceed with INNUca pipeline'
            runs['trueCoverage_ReMatCh'] = not_run
            runs['first_FastQC'] = not_run + ['NA']
            runs['Trimmomatic'] = not_run + ['NA']
            runs['second_Coverage'] = not_run
            runs['second_FastQC'] = not_run + ['NA']
            runs['Pear'] = not_run + ['NA']
            runs['SPAdes'] = not_run + ['NA']
            runs['Assembly_Mapping'] = not_run + ['NA']
            runs['Pilon'] = not_run
            runs['MLST'] = not_run + ['NA']

        if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
            if args.skipTrueCoverage or trueCoverage_config is None or (run_successfully_trueCoverage and pass_qc_trueCoverage):
                if args.skipFastQC or (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and runs['first_FastQC'][1])) is not False or args.fastQCproceed:
                    unassembled_pe_reads = None
                    assembled_se_reads = None
                    # Run Pear
                    if args.runPear:
                        print '--runPear set. Running Pear'
                        pearMinOverlap = pear.determine_minimum_overlap(args.pearMinOverlap, min_reads_length, max_reads_length)
                        run_successfully, pass_qc, time_taken, failing, unassembled_pe_reads, assembled_se_reads, pear_folder, warning = pear.runPear(fastq_files, threads, outdir, sampleName, fastq_encoding, trimmomatic_run_successfully, pearMinOverlap)
                        runs['Pear'] = [run_successfully, pass_qc, time_taken, failing, warning]
                    else:
                        runs['Pear'] = not_run + ['NA']

                    # Run SPAdes
                    if not args.skipSPAdes:
                        run_successfully, pass_qc, time_taken, failing, contigs_spades, warning = spades.runSpades(sampleName, outdir, threads, unassembled_pe_reads if unassembled_pe_reads is not None else fastq_files, args.spadesNotUseCareful, spadesMaxMemory, args.spadesMinCoverageAssembly, args.spadesMinContigsLength, genomeSize, args.spadesKmers, max_reads_length, args.spadesDefaultKmers, args.spadesMinKmerCovContigs, assembled_se_reads, args.saveExcludedContigs, args.maxNumberContigs)
                        runs['SPAdes'] = [run_successfully, pass_qc, time_taken, failing, warning]

                        if run_successfully:
                            contigs = contigs_spades

                            # Run Assembly Mapping check
                            bam_file = None
                            if not args.skipAssemblyMapping:
                                run_successfully, pass_qc, time_taken, failing, assembly_filtered, bam_file, assemblyMapping_folder, warning = assembly_mapping.runAssemblyMapping(fastq_files, contigs, threads, outdir, args.assemblyMinCoverageContigs, genomeSize, args.saveExcludedContigs, args.maxNumberContigs)
                                runs['Assembly_Mapping'] = [run_successfully, pass_qc, time_taken, failing, warning]

                                if run_successfully:
                                    contigs = assembly_filtered
                                    if not args.keepIntermediateAssemblies and os.path.isfile(contigs_spades) and contigs != contigs_spades:
                                        os.remove(contigs_spades)
                            else:
                                print '--skipAssemblyMapping set. Skipping Assembly Mapping check'
                                runs['Assembly_Mapping'] = skipped + ['NA']

                            # Run Pilon
                            if not args.skipPilon:
                                run_successfully, _, time_taken, failing, assembly_polished, pilon_folder = pilon.runPilon(jar_path_pilon, contigs, fastq_files, threads, outdir, jarMaxMemory, bam_file)
                                runs['Pilon'] = [run_successfully, None, time_taken, failing]

                                if run_successfully:
                                    contigs = assembly_polished
                                    if not args.keepIntermediateAssemblies and 'assembly_filtered' in locals() and os.path.isfile(assembly_filtered):
                                        os.remove(assembly_filtered)

                                if not args.pilonKeepFiles:
                                    utils.removeDirectory(pilon_folder)

                            else:
                                print '--skipPilon set. Skipping Pilon correction'
                                runs['Pilon'] = skipped

                            if 'assemblyMapping_folder' in locals():
                                utils.removeDirectory(assemblyMapping_folder)

                            print '\n' + 'Final assembly: ' + contigs
                            with open(os.path.join(outdir, 'final_assembly.txt'), 'wt') as writer:
                                writer.write(contigs + '\n')

                            # Run MLST
                            if not args.skipMLST:
                                run_successfully, pass_qc, time_taken, failing, warning = mlst.runMlst(contigs, scheme, outdir, species_genus, mlst_scheme_genus)
                                runs['MLST'] = [run_successfully, pass_qc, time_taken, failing, warning]
                            else:
                                print '--skipMLST set. Skipping MLST analysis'
                                runs['MLST'] = skipped + ['NA']
                        else:
                            print 'SPAdes did not run successfully! Skipping Pilon correction, Assembly Mapping check and MLST analysis'
                            runs['Assembly_Mapping'] = skipped + ['NA']
                            runs['Pilon'] = skipped
                            runs['MLST'] = skipped + ['NA']

                    else:
                        print '--skipSPAdes set. Skipping SPAdes, Pilon correction, Assembly Mapping check and MLST analysis'
                        runs['SPAdes'] = skipped + ['NA']
                        runs['Assembly_Mapping'] = skipped + ['NA']
                        runs['Pilon'] = skipped
                        runs['MLST'] = skipped + ['NA']
    else:
        print 'Moving to the next sample'
        for step in ('first_Coverage', 'trueCoverage_ReMatCh', 'first_FastQC', 'Trimmomatic', 'second_Coverage', 'second_FastQC', 'Pear', 'SPAdes', 'Assembly_Mapping', 'Pilon', 'MLST'):
            if step in ('Trimmomatic', 'first_FastQC', 'second_FastQC', 'Pear', 'SPAdes', 'Assembly_Mapping', 'MLST'):
                runs[step] = not_run + ['NA']
            else:
                runs[step] = not_run

    # Remove Pear directory
    if not args.pearKeepFiles and 'pear_folder' in locals():
        utils.removeDirectory(pear_folder)
    # Remove Trimmomatic directory with cleaned reads
    if not args.trimKeepFiles and 'trimmomatic_folder' in locals():
        utils.removeDirectory(trimmomatic_folder)

    # Check run
    run_successfully = all(runs[step][0] or runs[step][0] is None for step in runs)

    pass_fastqIntegrity = runs['FastQ_Integrity'][0]
    pass_cov = (runs['second_Coverage'][1] or (runs['second_Coverage'][1] is None and runs['first_Coverage'][1])) is not False
    pass_trueCov = runs['trueCoverage_ReMatCh'][1] is not False
    pass_fastqc = (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and runs['first_FastQC'][1])) is not False
    # pass_trimmomatic = runs['Trimmomatic'][1] is not False
    # pass_pear = runs['Pear'][1] is not False
    # pass_spades = runs['SPAdes'][1] is not False or runs['Assembly_Mapping'][1] is True
    pass_spades = runs['SPAdes'][1] is not False
    pass_assemblyMapping = runs['Assembly_Mapping'][1] is not False
    pass_pilon = runs['Pilon'][0] is not False
    pass_mlst = runs['MLST'][1] is not False
    pass_qc = all([pass_fastqIntegrity, pass_cov, pass_trueCov, pass_fastqc, pass_spades, pass_assemblyMapping, pass_pilon, pass_mlst])

    return run_successfully, pass_qc, runs


if __name__ == "__main__":
    main()
