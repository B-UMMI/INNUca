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
import modules.fastQintegrity as fastQintegrity
import modules.estimated_coverage as coverage
import modules.fastqc as fastqc
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
	outdir = os.path.abspath(os.path.join(args.outdir, ''))
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
	programs_version_dictionary = {}
	if not args.skipEstimatedCoverage:
		programs_version_dictionary['bunzip2'] = ['--version', '>=', '1.0.6']
		programs_version_dictionary['gunzip'] = ['--version', '>=', '1.6']
	if not (args.skipFastQC and args.skipTrimmomatic):
		programs_version_dictionary['java'] = ['-version', '>=', '1.8']
	if not args.skipFastQC:
		programs_version_dictionary['fastqc'] = ['--version', '==', '0.11.5']
	if not args.skipTrimmomatic:
		programs_version_dictionary['trimmomatic-0.36.jar'] = ['-version', '==', '0.36']
	if not args.skipSPAdes:
		programs_version_dictionary['spades.py'] = ['--version', '>=', '3.7.1']
	if not args.skipMLST and not args.skipSPAdes:
		programs_version_dictionary['mlst'] = ['--version', '>=', '2.4']
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

	# runningTime in seconds
	# fileSize in bytes

	number_samples_successfully = 0
	number_samples_pass = 0

	steps = ['FastQ_Integrity', 'first_Coverage', 'first_FastQC', 'second_Coverage', 'second_FastQC', 'SPAdes', 'MLST']

	sample_lines = []
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

			failures = []

			for step in steps:
				max_fails = len(list(run_report[step][3].values()))

				if fail_reasons.count(False) < len(max_fails):

					failures.append('#' + step)

					for key, fail_reasons in run_report[step][3].items():
						if isinstance(fail_reasons, bool) and not fail_reasons:
							continue
						else:
							failures.append('>' + str(key))

							if isinstance(fail_reasons, (list, tuple)):
								for reasons in fail_reasons:
									failures.append(str(reasons))
							else:
								failures.append(str(fail_reasons))

			writer_failReport.write( '\n'.join(failures) )

		# Save runs statistics
		if run_successfully:
			number_samples_successfully += 1
		if pass_qc:
			number_samples_pass += 1

		# Get raw reads files size
		fileSize = sum(os.path.getsize(fastq) for fastq in fastq_files)

		# Remove sample directory if it was created during the process
		if removeCreatedSamplesDirectories:
			utils.removeDirectory(os.path.join(inputDirectory, sample, ''))

		print 'END ' + sample + ' analysis'
		time_taken = utils.runTime(sample_start_time)

		sample_lines.append(utils.sampleReportLine(run_report))

	# Prepare run report file

	samples_report_path = os.path.join(outdir, 'samples_report.' + time_str + '.tab')

	with open(samples_report_path, 'w') as sample_report:

		out = csv.writer(sample_report, delimiter='\t')

		header = utils.build_header(steps)

		out.writerow(header)

		for line in sample_lines:
			out.writerow(line)
	# Run report
	print '\n' + 'END INNUca.py'
	print '\n' + str(number_samples_successfully) + ' samples out of ' + str(len(samples)) + ' run successfully'
	print '\n' + str(number_samples_pass) + ' samples out of ' + str(number_samples_successfully) + ' (run successfully) PASS INNUca.py analysis'
	time_taken = utils.runTime(general_start_time)

	# Check whether INNUca.py run at least one sample successfully
	if number_samples_successfully == 0:
		sys.exit('No samples run successfully!')


def run_INNUca(sampleName, outdir, fastq_files, args, script_path):
	threads = args.threads
	adaptersFasta = args.adapters
	if adaptersFasta is not None:
		adaptersFasta = os.path.abspath(adaptersFasta.name)
	genomeSize = args.genomeSizeExpectedMb
	maximumReadsLength = None

	runs = {}

	# Run FastQ integrity check
	program_start_time = time.time()
	print 'RUNNING FastQ integrity check'
	not_corruption_found, failing = fastQintegrity.runFastQintegrity(fastq_files, threads, outdir)
	print 'END FastQ integrity check'
	time_taken = utils.runTime(program_start_time)
	runs['FastQ_Integrity'] = [not_corruption_found, None, time_taken, failing]

	if not_corruption_found:
		# Run first Estimated Coverage
		if not args.skipEstimatedCoverage:

			program_start_time = time.time()
			print 'RUNNING First Estimated Coverage analysis'

			# Check whether the Estimated Coverage output is already present
			report_file = os.path.join(outdir, 'coverage_report.txt')

			if os.path.isfile(report_file):
				os.remove(report_file)

			# Run getEstimatedCoverage
			run_successfully, pass_qc, failing = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir)
			print 'END First Estimated Coverage analysis'

			time_taken = utils.runTime(program_start_time)

			runs['first_Coverage'] = [run_successfully, pass_qc, time_taken, failing]
		else:
			print '--skipEstimatedCoverage set. Skipping First Estimated Coverage analysis'
			runs['first_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]

		# Run first FastQC
		nts2clip_based_ntsContent = None

		if not args.skipFastQC:

			program_start_time = time.time()
			print 'RUNNING First FastQC analysis'

			run_successfully, pass_qc, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
			print 'END First FastQC analysis'
			time_taken = utils.runTime(program_start_time)

			runs['first_FastQC'] = [run_successfully, pass_qc, time_taken, failing]

		else:
			print '--skipFastQC set. Skipping First FastQC analysis'
			runs['first_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]

		# Run Trimmomatic
		if not args.skipTrimmomatic:
			program_start_time = time.time()
			print 'RUNNING Trimmomatic'

			run_successfully, paired_reads, trimmomatic_folder, failing = trimmomatic.runTrimmomatic(sampleName, outdir, threads, adaptersFasta, script_path, args.doNotSearchAdapters, fastq_files, maximumReadsLength, args.doNotTrimCrops, args.trimCrop, args.trimHeadCrop, args.trimLeading, args.trimTrailing, args.trimSlidingWindow, args.trimMinLength, nts2clip_based_ntsContent)
			print 'END Trimmomatic'

			time_taken = utils.runTime(program_start_time)

			runs['Trimmomatic'] = [run_successfully, None, time_taken, failing]

			if run_successfully:
				fastq_files = paired_reads

				# Run second Estimated Coverage
				if not args.skipEstimatedCoverage:
					program_start_time = time.time()
					print 'RUNNING Second Estimated Coverage analysis'

					run_successfully, pass_qc, failing = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir)
					print 'END Second Estimated Coverage analysis'
					time_taken = utils.runTime(program_start_time)

					runs['second_Coverage'] = [run_successfully, pass_qc, time_taken, failing]
				else:
					print '--skipEstimatedCoverage set. Skipping Second Estimated Coverage analysis'
					runs['second_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]

				# Run second FastQC
				if not args.skipFastQC:
					program_start_time = time.time()
					print 'RUNNING Second FastQC analysis'

					run_successfully, pass_qc, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
					print 'END Second FastQC analysis'
					time_taken = utils.runTime(program_start_time)

					runs['second_FastQC'] = [run_successfully, pass_qc, time_taken, failing]
				else:
					print '--skipFastQC set. Skipping Second FastQC analysis'
					runs['second_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]
			else:
				print 'Trimmomatic did not run successfully! Skipping Second Estimated Coverage analysis and FastQC analysis'
				runs['second_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]
				runs['second_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]

		else:
			print '--skipTrimmomatic set. Skipping Trimmomatic, but also Second FastQC analysis and Second Estimated Coverage analysis'
			runs['Trimmomatic'] = [None, None, 0, {'sample': 'Skipped'}]
			runs['second_Coverage'] = [None, None, 0, {'sample': 'Skipped'}]
			runs['second_FastQC'] = [None, None, 0, {'sample': 'Skipped'}]

		# Run SPAdes
		if not args.skipSPAdes:
			program_start_time = time.time()
			print 'RUNNING SPAdes'

			run_successfully, pass_qc, failing, contigs = spades.runSpades(sampleName, outdir, threads, fastq_files, args.spadesNotUseCareful, args.spadesMaxMemory[0], args.spadesMinCoverage[0], args.spadesMinContigsLength[0], genomeSize, args.spadesKmers, maximumReadsLength, args.spadesSaveReport, args.spadesDefaultKmers)
			print 'END SPAdes'
			time_taken = utils.runTime(program_start_time)

			runs['SPAdes'] = [run_successfully, pass_qc, time_taken, failing]

			if run_successfully:
				# Run MLST
				if not args.skipMLST:
					program_start_time = time.time()
					print 'RUNNING MLST analysis'

					run_successfully, pass_qc, failing = mlst.runMlst(contigs, args.speciesExpected[0], outdir)
					print 'END MLST analysis'
					time_taken = utils.runTime(program_start_time)

					runs['MLST'] = [run_successfully, pass_qc, time_taken, failing]
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
	else:
		print 'Moving to the next sample'
		runs['first_Coverage'] = [None, None, 0, {'sample': 'Not run'}]
		runs['first_FastQC'] = [None, None, 0, {'sample': 'Not run'}]
		runs['Trimmomatic'] = [None, None, 0, {'sample': 'Not run'}]
		runs['second_Coverage'] = [None, None, 0, {'sample': 'Not run'}]
		runs['second_FastQC'] = [None, None, 0, {'sample': 'Not run'}]
		runs['SPAdes'] = [None, None, 0, {'sample': 'Not run'}]
		runs['MLST'] = [None, None, 0, {'sample': 'Not run'}]

	# Remove Trimmomatic directory with cleaned reads
	if not args.trimKeepFiles:
		try:
			utils.removeDirectory(trimmomatic_folder)
		except:
			print 'It is not possible to remove Trimmomatic directory because Trimmomatic did not run'

	# Check run
	run_successfully = True
	pass_qc = True

	for step in runs:
		if step != 'FastQ_Integrity' and runs[step][0] is False:
			run_successfully = False
			break

	if runs['FastQ_Integrity'][0] is False:
		pass_qc = False

	if runs['second_Coverage'][1] is False:
		pass_qc = False
	elif runs['second_Coverage'][1] is None and runs['first_Coverage'][1] is False:
		pass_qc = False

	if runs['second_FastQC'][1] is False:
		pass_qc = False
	elif runs['second_FastQC'][1] is None and runs['first_FastQC'][1] is False:
		pass_qc = False

	if runs['SPAdes'][1] is False:
		pass_qc = False

	if runs['MLST'][1] is False:
		pass_qc = False

	return run_successfully, pass_qc, runs


if __name__ == "__main__":
	main()
