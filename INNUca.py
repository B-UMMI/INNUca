#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
INNUca - Reads Control and Assembly
INNUca.py - INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search
<https://github.com/B-UMMI/INNUca>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: July 01, 2016

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
import modules.fastqc as fastqc
import modules.estimated_coverage as coverage
import modules.trimmomatic as trimmomatic
import modules.spades as spades
import modules.mlst as mlst
import time
import os
import sys


def main():
	version = '1.4'
	args = utils.parseArguments(version)

	general_start_time = time.time()
	time_str = time.strftime("%Y%m%d-%H%M%S")

	# Check if output directory exists
	outdir = os.path.abspath(os.path.join(args.outdir[0], ''))
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	# Start logger
	sys.stdout = utils.Logger(outdir, time_str)

	print '\n' + '==========> INNUca.py <=========='
	print '\n' + 'Program start: ' + time.ctime()

	# Tells where the logfile will be stored
	print '\n' + 'LOGFILE:'
	print sys.stdout.getLogFile()

	# Print command
	print '\n' + 'COMMAND:'
	script_path = os.path.abspath(sys.argv[0])
	print sys.executable + ' ' + script_path + ' ' + ' '.join(sys.argv[1:])

	# Print directory where programme was lunch
	print '\n' + 'PRESENT DIRECTORY :'
	print os.getcwd()

	# Print program version
	print '\n' + 'VERSION INNUca.py:'
	utils.scriptVersionGit(version, os.getcwd(), script_path)

	# Set and print PATH variable
	utils.setPATHvariable(args.doNotUseProvidedSoftware, script_path)

	# Check programms
	programs_version_dictionary = {'java': ['-version', '>=', '1.8'], 'bunzip2': ['--version', '>=', '1.0.6'], 'gunzip': ['--version', '>=', '1.6'], 'fastqc': ['--version', '==', '0.11.5'], 'trimmomatic-0.36.jar': ['-version', '==', '0.36'], 'spades.py': ['--version', '>=', '3.7.1'], 'mlst': ['--version', '>=', '2.4']}
	missingPrograms = utils.checkPrograms(programs_version_dictionary)
	if len(missingPrograms) > 0:
		sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))

	# Check if input directory exists with fastq files and store samples name that have fastq files
	inputDirectory = os.path.abspath(os.path.join(args.inputDirectory[0], ''))
	pairEnd_filesSeparation_list = args.pairEnd_filesSeparation
	print ''
	samples, removeCreatedSamplesDirectories = utils.checkSetInputDirectory(inputDirectory, outdir, pairEnd_filesSeparation_list)

	# Start running the analysis
	print '\n' + 'RUNNING INNUca.py'

	steps = ['first_Coverage', 'first_FastQC', 'Trimmomatic', 'second_Coverage', 'second_FastQC', 'SPAdes', 'MLST']

	# Prepare run report file
	samples_report = open(os.path.join(outdir, str('samples_report.' + time_str + '.tab')), 'wt')
	# runningTime in seconds
	# fileSize in bytes
	samples_report.write('#samples' + '\t' + 'samples_runSuccessfully' + '\t' + 'samples_passQC' + '\t' + 'samples_runningTime' + '\t' + 'samples_fileSize' + '\t')
	for step in steps:
		if step == 'Trimmomatic':
			samples_report.write(str(step + '_runSuccessfully') + '\t' + str(step + '_runningTime') + '\t')
		else:
			samples_report.write(str(step + '_runSuccessfully') + '\t' + str(step + '_passQC') + '\t' + str(step + '_runningTime'))
			if step != steps[len(steps) - 1]:
				samples_report.write('\t')
			else:
				samples_report.write('\n')
	samples_report.flush()

	number_samples_successfully = 0
	number_samples_pass = 0

	# Run comparisons for each sample
	for sample in samples:
		sample_start_time = time.time()

		print '\n' + 'Sample: ' + sample

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

		print 'The following files will be used:'
		print str(fastq_files)

		# Run INNUca.py analysis
		run_successfully, pass_qc, run_report = run_INNUca(sample, sample_outdir, fastq_files, args, script_path)

		# Save sample fail report
		with open(os.path.join(sample_outdir, 'fail_report.txt'), 'wt') as writer_failReport:
			for step in steps:
				fail_reasons = list(run_report[step][3].values())
				if fail_reasons.count(False) == len(fail_reasons):
					continue
				else:
					writer_failReport.write('#' + step + '\n')
					for key, fail_reasons in run_report[step][3].items():
						if isinstance(fail_reasons, bool) and not fail_reasons:
							continue
						else:
							writer_failReport.write('>' + str(key) + '\n')
							if isinstance(fail_reasons, (list, tuple)):
								for reasons in fail_reasons:
									writer_failReport.write(str(reasons) + '\n')
									writer_failReport.flush()
							else:
								writer_failReport.write(str(fail_reasons) + '\n')
								writer_failReport.flush()

		# Save runs statistics
		if run_successfully:
			number_samples_successfully = number_samples_successfully + 1
		if pass_qc:
			number_samples_pass = number_samples_pass + 1

		# Get raw reads files size
		fileSize = 0
		for fastq in fastq_files:
			fileSize = fileSize + os.path.getsize(fastq)

		# Remove sample directory if it was created during the process
		if removeCreatedSamplesDirectories:
			utils.removeDirectory(os.path.join(inputDirectory, sample, ''))

		print 'END ' + sample + ' analysis'
		time_taken = utils.runTime(sample_start_time)

		# Save run report
		samples_report.write(sample + '\t' + str(run_successfully) + '\t' + ('PASS' if pass_qc else 'FAIL') + '\t' + str(time_taken) + '\t' + str(fileSize) + '\t')
		for step in steps:
			if step == 'Trimmomatic':
				samples_report.write(str(run_report[step][0]) + '\t' + str(run_report[step][2]) + '\t')
			else:
				samples_report.write(str(run_report[step][0]) + '\t' + ('PASS' if run_report[step][1] else 'FAIL') + '\t' + str(run_report[step][2]))
				if step != steps[len(steps) - 1]:
					samples_report.write('\t')
				else:
					samples_report.write('\n')
		samples_report.flush()

	samples_report.close()

	# Run report
	print '\n' + 'END INNUca.py'
	print '\n' + str(number_samples_successfully) + ' samples out of ' + str(len(samples)) + ' run successfully'
	print '\n' + str(number_samples_pass) + ' samples out of ' + str(number_samples_successfully) + ' (run successfully) PASS INNUca.py analysis'
	time_taken = utils.runTime(general_start_time)

	# Check whether INNUca.py run at least one sample successfully
	if number_samples_successfully == 0:
		sys.exit('No samples run successfully!')


def run_INNUca(sampleName, outdir, fastq_files, args, script_path):
	threads = args.threads[0]
	adaptersFasta = args.adapters[0]
	if adaptersFasta is not None:
		adaptersFasta = os.path.abspath(adaptersFasta.name)
	genomeSize = args.genomeSizeExpectedMb[0]
	maximumReadsLength = None

	runs = {}

	# Run first Estimated Coverage
	if not args.skipEstimatedCoverage:
		program_start_time = time.time()
		print 'RUNNING First Estimated Coverage analysis'
		program = []
		# Check whether the Estimated Coverage output is already present
		report_file = os.path.join(outdir, 'coverage_report.txt')
		if os.path.isfile(report_file):
			os.remove(report_file)
		# Run getEstimatedCoverage
		run_successfully, pass_qc, failing = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir)
		print 'END First Estimated Coverage analysis'
		time_taken = utils.runTime(program_start_time)
		program.append(run_successfully)
		program.append(pass_qc)
		program.append(time_taken)
		program.append(failing)
		runs['first_Coverage'] = []
		runs['first_Coverage'] = program
	else:
		print '--skipEstimatedCoverage set. Skipping First Estimated Coverage analysis'
		runs['first_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]

	# Run first FastQC
	nts2clip_based_ntsContent = None
	if not args.skipFastQC:
		program_start_time = time.time()
		print 'RUNNING First FastQC analysis'
		program = []
		run_successfully, pass_qc, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
		print 'END First FastQC analysis'
		time_taken = utils.runTime(program_start_time)
		program.append(run_successfully)
		program.append(pass_qc)
		program.append(time_taken)
		program.append(failing)
		runs['first_FastQC'] = []
		runs['first_FastQC'] = program
	else:
		print '--skipFastQC set. Skipping First FastQC analysis'
		runs['first_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]

	# Run Trimmomatic
	if not args.skipTrimmomatic:
		program_start_time = time.time()
		print 'RUNNING Trimmomatic'
		program = []
		run_successfully, paired_reads, trimmomatic_folder, failing = trimmomatic.runTrimmomatic(sampleName, outdir, threads, adaptersFasta, script_path, args.doNotSearchAdapters, fastq_files, maximumReadsLength, args.doNotTrimCrops, args.trimCrop, args.trimHeadCrop, args.trimLeading[0], args.trimTrailing[0], args.trimSlidingWindow[0], args.trimMinLength[0], nts2clip_based_ntsContent)
		print 'END Trimmomatic'
		time_taken = utils.runTime(program_start_time)
		program.append(run_successfully)
		program.append(None)
		program.append(time_taken)
		program.append(failing)
		runs['Trimmomatic'] = []
		runs['Trimmomatic'] = program

		if run_successfully:
			fastq_files = paired_reads

			# Run second Estimated Coverage
			if not args.skipEstimatedCoverage:
				program_start_time = time.time()
				print 'RUNNING Second Estimated Coverage analysis'
				program = []
				run_successfully, pass_qc, failing = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir)
				print 'END Second Estimated Coverage analysis'
				time_taken = utils.runTime(program_start_time)
				program.append(run_successfully)
				program.append(pass_qc)
				program.append(time_taken)
				program.append(failing)
				runs['second_Coverage'] = []
				runs['second_Coverage'] = program
			else:
				print '--skipEstimatedCoverage set. Skipping Second Estimated Coverage analysis'
				runs['second_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]

			# Run second FastQC
			if not args.skipFastQC:
				program_start_time = time.time()
				print 'RUNNING Second FastQC analysis'
				program = []
				run_successfully, pass_qc, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
				print 'END Second FastQC analysis'
				time_taken = utils.runTime(program_start_time)
				program.append(run_successfully)
				program.append(pass_qc)
				program.append(time_taken)
				program.append(failing)
				runs['second_FastQC'] = []
				runs['second_FastQC'] = program
			else:
				print '--skipFastQC set. Skipping Second FastQC analysis'
				runs['second_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]
		else:
			print 'Trimmomatic did not run successfully! Skipping Second Estimated Coverage analysis and FastQC analysis'
			runs['second_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]
			runs['second_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]

	else:
		print '--skipTrimmomatic set. Skipping Trimmomatic'
		runs['Trimmomatic'] = [None, None, 0, {'sample': 'Skipped'}]

	# Run SPAdes
	if not args.skipSPAdes:
		program_start_time = time.time()
		print 'RUNNING SPAdes'
		program = []
		run_successfully, pass_qc, failing, contigs = spades.runSpades(sampleName, outdir, threads, fastq_files, args.spadesNotUseCareful, args.spadesMaxMemory[0], args.spadesMinCoverage[0], args.spadesMinContigsLength[0], genomeSize, args.spadesKmers, maximumReadsLength, args.spadesSaveReport)
		print 'END SPAdes'
		time_taken = utils.runTime(program_start_time)
		program.append(run_successfully)
		program.append(pass_qc)
		program.append(time_taken)
		program.append(failing)
		runs['SPAdes'] = []
		runs['SPAdes'] = program

		if run_successfully:
			# Run MLST
			if not args.skipMLST:
				program_start_time = time.time()
				print 'RUNNING MLST analysis'
				program = []
				run_successfully, pass_qc, failing = mlst.runMlst(contigs, args.speciesExpected[0], outdir)
				print 'END MLST analysis'
				time_taken = utils.runTime(program_start_time)
				program.append(run_successfully)
				program.append(pass_qc)
				program.append(time_taken)
				program.append(failing)
				runs['MLST'] = []
				runs['MLST'] = program
			else:
				print '--skipMLST set. Skipping MLST analysis'
				runs['MLST'] = [None, None, 0, {'sample': 'Skipped'}]
		else:
			print 'SPAdes did not run successfully! Skipping MLST analysis'
			runs['MLST'] = [None, None, 0, {'sample': 'Skipped'}]

	else:
		print '--skipSPAdes set. Skipping SPAdes and MLST analysis'
		runs['SPAdes'] = [None, None, 0, {'sample': 'Skipped'}]
		runs['MLST'] = [None, None, 0, {'sample': 'Skipped'}]

	# Remove Trimmomatic directory with cleaned reads
	if not args.trimKeepFiles:
		try:
			utils.removeDirectory(trimmomatic_folder)
		except:
			print 'It is not possible to remove Trimmomatic directory ' + trimmomatic_folder + ' because Trimmomatic did not run'

	# Check run
	run_successfully = True
	pass_qc = True

	for step in runs:
		if not runs[step][0]:
			run_successfully = False

	if not runs['second_Coverage'][1]:
		pass_qc = False
	elif runs['second_Coverage'][1] is None and not runs['first_Coverage'][1]:
		pass_qc = False
	elif runs['second_Coverage'][1] is None and not runs['first_Coverage'][1]:
		pass_qc = False

	if not runs['second_FastQC'][1]:
		pass_qc = False
	elif runs['second_FastQC'][1] is None and not runs['first_FastQC'][1]:
		pass_qc = False
	elif runs['second_FastQC'][1] is None and not runs['first_FastQC'][1]:
		pass_qc = False

	if not runs['SPAdes'][1] or runs['SPAdes'][1] is None:
		pass_qc = False

	if not runs['MLST'][1] or runs['MLST'][1] is None:
		pass_qc = False

	return run_successfully, pass_qc, runs


if __name__ == "__main__":
	main()
