#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
INNUca - Reads Control and Assembly
INNUca.py - INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search
<https://github.com/B-UMMI/INNUca>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: December 23, 2016

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
import modules.pilon as pilon
import modules.mlst as mlst
import modules.assembly_mapping as assembly_mapping
import modules.trueCoverage_rematch as trueCoverage
import time
import os
import sys


def main():
	version = '2.1'
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

	# Get CPU information
	utils.get_cpu_information(outdir, time_str)

	# Set and print PATH variable
	utils.setPATHvariable(args.doNotUseProvidedSoftware, script_path)

	# Check programms
	programs_version_dictionary = {}
	programs_version_dictionary['gunzip'] = ['--version', '>=', '1.6']
	if (not args.skipTrueCoverage or (not args.skipPilon and not args.skipSPAdes)):
		programs_version_dictionary['bowtie2'] = ['--version', '>=', '2.2.9']
		programs_version_dictionary['samtools'] = ['--version', '==', '1.3.1']
	if not (args.skipFastQC and args.skipTrimmomatic and (args.skipPilon or args.skipSPAdes)):
		# programs_version_dictionary['java'] = ['-version', '>=', '1.8']
		programs_version_dictionary['java'] = [None, '>=', '1.8']  # For OpenJDK compatibility
	if not args.skipFastQC:
		programs_version_dictionary['fastqc'] = ['--version', '==', '0.11.5']
	if not args.skipTrimmomatic:
		programs_version_dictionary['trimmomatic-0.36.jar'] = ['-version', '==', '0.36']
	if not args.skipPear or not args.skipSPAdes:
		programs_version_dictionary['pear'] = ['--version', '>=', '0.9.10']
	if not args.skipSPAdes:
		programs_version_dictionary['spades.py'] = ['--version', '>=', '3.9.0']
	if not args.skipPilon and not args.skipSPAdes:
		programs_version_dictionary['pilon-1.18.jar'] = ['--version', '==', '1.18']
	if not args.skipMLST and not args.skipSPAdes:
		programs_version_dictionary['mlst'] = ['--version', '>=', '2.4']
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

	# Check if input directory exists with fastq files and store samples name that have fastq files
	inputDirectory = os.path.abspath(os.path.join(args.inputDirectory, ''))
	# pairEnd_filesSeparation_list = args.pairEnd_filesSeparation
	pairEnd_filesSeparation_list = None
	print ''
	samples, removeCreatedSamplesDirectories, indir_same_outdir = utils.checkSetInputDirectory(inputDirectory, outdir, pairEnd_filesSeparation_list)

	# Start running the analysis
	print '\n' + 'RUNNING INNUca.py'

	# Prepare run report file
	samples_report_path = os.path.join(outdir, 'samples_report.' + time_str + '.tab')
	utils.start_sample_report_file(samples_report_path)

	number_samples_successfully = 0
	number_samples_pass = 0

	# Get MLST scheme to use
	scheme = 'unknown'
	if not args.skipMLST and not args.skipSPAdes:
		scheme = mlst.getScheme(args.speciesExpected)

	# Get path to blastn
	mlst.getBlastPath()

	# Get trueCoverage_ReMatCh settings
	trueCoverage_config = None
	if not args.skipTrueCoverage:
		trueCoverage_reference = None
		trueCoverage_config_file = None
		trueCoverage_config = None

		if args.trueConfigFile is None:
			print 'No trueCoverage_ReMatCh config file was provided. Search for default files'
			trueCoverage_config_file, trueCoverage_reference = trueCoverage.check_existing_default_config(args.speciesExpected, script_path)
		else:
			trueCoverage_config_file = args.trueConfigFile.name

		if trueCoverage_config_file is not None:
			trueCoverage_config = trueCoverage.parse_config(trueCoverage_config_file)
		if args.trueConfigFile is None and trueCoverage_config is not None:
			trueCoverage_config['reference_file'] = trueCoverage_reference

		if trueCoverage_config is not None:
			print 'The following trueCoverage_ReMatCh config file will be used: ' + trueCoverage_config_file
			print 'The following trueCoverage_ReMatCh reference file will be used: ' + trueCoverage_config['reference_file'] + '\n'
		else:
			print 'No trueCoverage_ReMatCh config file was found'

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

		print 'The following files will be used:'
		print str(fastq_files) + '\n'

		# Run INNUca.py analysis
		run_successfully, pass_qc, run_report = run_INNUca(sample, sample_outdir, fastq_files, args, script_path, scheme, spadesMaxMemory, jar_path_trimmomatic, jar_path_pilon, jarMaxMemory, trueCoverage_config)

		# Save sample fail report
		fail_report_path = os.path.join(sample_outdir, 'fail_report.txt')
		utils.write_fail_report(fail_report_path, run_report)

		# Save runs statistics
		if run_successfully:
			number_samples_successfully += 1
		if pass_qc:
			number_samples_pass += 1

		# Get raw reads files size
		fileSize = sum(os.path.getsize(fastq) for fastq in fastq_files)

		# Remove sample directory if it was created during the process
		if removeCreatedSamplesDirectories and not indir_same_outdir:
			utils.removeDirectory(os.path.join(inputDirectory, sample, ''))

		print 'END ' + sample + ' analysis'
		time_taken = utils.runTime(sample_start_time)

		# Save run report
		utils.write_sample_report(samples_report_path, sample, run_successfully, pass_qc, time_taken, fileSize, run_report)

	# Run report
	print '\n' + 'END INNUca.py'
	print '\n' + str(number_samples_successfully) + ' samples out of ' + str(len(samples)) + ' run successfully'
	print '\n' + str(number_samples_pass) + ' samples out of ' + str(number_samples_successfully) + ' (run successfully) PASS INNUca.py analysis'
	time_taken = utils.runTime(general_start_time)
	del time_taken

	# Check whether INNUca.py run at least one sample successfully
	if number_samples_successfully == 0:
		sys.exit('No samples run successfully!')


def run_INNUca(sampleName, outdir, fastq_files, args, script_path, scheme, spadesMaxMemory, jar_path_trimmomatic, jar_path_pilon, jarMaxMemory, trueCoverage_config):
	threads = args.threads
	adaptersFasta = args.adapters
	if adaptersFasta is not None:
		adaptersFasta = os.path.abspath(adaptersFasta.name)
	genomeSize = args.genomeSizeExpectedMb
	maximumReadsLength = None
	skipped = [None, None, 0, {'sample': 'Skipped'}]
	not_run = [None, None, 0, {'sample': 'Not run'}]

	runs = {}

	# Run FastQ integrity check
	not_corruption_found, _, time_taken, failing, fastq_encoding = fastQintegrity.runFastQintegrity(fastq_files, threads, outdir)
	runs['FastQ_Integrity'] = [not_corruption_found, None, time_taken, failing]

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
			run_successfully_estimatedCoverage, pass_qc, time_taken, failing, estimatedCoverage = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir, threads)
			runs['first_Coverage'] = [run_successfully_estimatedCoverage, pass_qc, time_taken, failing]
		else:
			print '--skipEstimatedCoverage set. Skipping First Estimated Coverage analysis'
			runs['first_Coverage'] = skipped

		trimmomatic_run_successfully = False

		if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
			if not args.skipTrueCoverage and trueCoverage_config is not None:
				# Run True Coverage
				run_successfully_trueCoverage, pass_qc_trueCoverage, time_taken, failing = trueCoverage.runTrueCoverage(fastq_files, trueCoverage_config['reference_file'], threads, outdir, trueCoverage_config['length_extra_seq'], trueCoverage_config['minimum_depth_presence'], trueCoverage_config['minimum_depth_call'], trueCoverage_config['minimum_depth_frequency_dominant_allele'], trueCoverage_config['minimum_gene_coverage'], trueCoverage_config['maximum_number_absent_genes'], trueCoverage_config['maximum_number_genes_multiple_alleles'], trueCoverage_config['minimum_read_coverage'])
				runs['trueCoverage_ReMatCh'] = [run_successfully_trueCoverage, pass_qc_trueCoverage, time_taken, failing]
			else:
				print '\n' + '--skipTrueCoverage set. Skipping True coverage analysis'
				runs['trueCoverage_ReMatCh'] = skipped

			if args.skipTrueCoverage or trueCoverage_config is None or (run_successfully_trueCoverage and pass_qc_trueCoverage):
				# Run first FastQC
				nts2clip_based_ntsContent = None
				if not args.skipFastQC:
					run_successfully, pass_qc, time_taken, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
					runs['first_FastQC'] = [run_successfully, pass_qc, time_taken, failing]
				else:
					print '--skipFastQC set. Skipping First FastQC analysis'
					runs['first_FastQC'] = skipped

				# Run Trimmomatic
				if not args.skipTrimmomatic:
					run_successfully, not_empty_fastq, time_taken, failing, paired_reads, trimmomatic_folder, fileSize = trimmomatic.runTrimmomatic(jar_path_trimmomatic, sampleName, outdir, threads, adaptersFasta, script_path, args.doNotSearchAdapters, fastq_files, maximumReadsLength, args.doNotTrimCrops, args.trimCrop, args.trimHeadCrop, args.trimLeading, args.trimTrailing, args.trimSlidingWindow, args.trimMinLength, nts2clip_based_ntsContent, jarMaxMemory, fastq_encoding)
					runs['Trimmomatic'] = [run_successfully, not_empty_fastq, time_taken, failing, fileSize]
					trimmomatic_run_successfully = run_successfully

					if run_successfully and not_empty_fastq:
						fastq_files = paired_reads

						# Run second Estimated Coverage
						if not args.skipEstimatedCoverage:
							run_successfully_estimatedCoverage, pass_qc, time_taken, failing, estimatedCoverage = coverage.getEstimatedCoverage(fastq_files, genomeSize, outdir, threads)
							runs['second_Coverage'] = [run_successfully_estimatedCoverage, pass_qc, time_taken, failing]
						else:
							print '--skipEstimatedCoverage set. Skipping Second Estimated Coverage analysis'
							runs['second_Coverage'] = skipped

						if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
							# Run second FastQC
							if not args.skipFastQC:
								run_successfully, pass_qc, time_taken, failing, maximumReadsLength, nts2clip_based_ntsContent = fastqc.runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files)
								runs['second_FastQC'] = [run_successfully, pass_qc, time_taken, failing]
							else:
								print '--skipFastQC set. Skipping Second FastQC analysis'
								runs['second_FastQC'] = skipped
						else:
							print '\n' + 'Estimated coverage is too lower (< ' + str(args.estimatedMinimumCoverage) + 'x). This sample will not proceed with INNUca pipeline'
							runs['second_FastQC'] = not_run
							runs['SPAdes'] = not_run
							runs['Pilon'] = not_run
							runs['Assembly_Mapping'] = not_run
							runs['MLST'] = not_run
					else:
						print 'Trimmomatic did not run successfully or return zero reads! Skipping Second Estimated Coverage analysis and FastQC analysis'
						runs['second_Coverage'] = skipped
						runs['second_FastQC'] = skipped

				else:
					print '--skipTrimmomatic set. Skipping Trimmomatic, but also Second FastQC analysis and Second Estimated Coverage analysis'
					runs['Trimmomatic'] = skipped + ['NA']
					runs['second_Coverage'] = skipped
					runs['second_FastQC'] = skipped
			else:
				print '\n' + 'This sample does not pass True Coverage module QA/QC. This sample will not proceed with INNUca pipeline'
				runs['first_FastQC'] = not_run
				runs['Trimmomatic'] = not_run + ['NA']
				runs['second_Coverage'] = not_run
				runs['second_FastQC'] = not_run
				runs['Pear'] = not_run
				runs['SPAdes'] = not_run
				runs['Pilon'] = not_run
				runs['Assembly_Mapping'] = not_run
				runs['MLST'] = not_run

		else:
			print '\n' + 'Estimated coverage is too lower (< ' + str(args.estimatedMinimumCoverage) + 'x). This sample will not proceed with INNUca pipeline'
			runs['trueCoverage_ReMatCh'] = not_run
			runs['first_FastQC'] = not_run
			runs['Trimmomatic'] = not_run + ['NA']
			runs['second_Coverage'] = not_run
			runs['second_FastQC'] = not_run
			runs['Pear'] = not_run
			runs['SPAdes'] = not_run
			runs['Pilon'] = not_run
			runs['Assembly_Mapping'] = not_run
			runs['MLST'] = not_run

		if args.skipEstimatedCoverage or (run_successfully_estimatedCoverage and not estimatedCoverage < args.estimatedMinimumCoverage):
			if args.skipTrueCoverage or trueCoverage_config is None or (run_successfully_trueCoverage and pass_qc_trueCoverage):
				unassembled_pe_reads = None
				assembled_se_reads = None
				# Run Pear
				if not args.skipPear or not args.skipSPAdes:
					run_successfully, pass_qc, time_taken, failing, unassembled_pe_reads, assembled_se_reads, pear_folder = pear.runPear(fastq_files, threads, outdir, sampleName, fastq_encoding, trimmomatic_run_successfully)
					runs['Pear'] = [run_successfully, pass_qc, time_taken, failing]
				else:
					print '--skipPear or --skipSPAdes set. Skipping Pear'
					runs['Pear'] = skipped

				# Run SPAdes
				if not args.skipSPAdes:
					run_successfully, pass_qc, time_taken, failing, contigs_spades = spades.runSpades(sampleName, outdir, threads, unassembled_pe_reads if unassembled_pe_reads is not None else fastq_files, args.spadesNotUseCareful, spadesMaxMemory, args.spadesMinCoverageAssembly, args.spadesMinContigsLength, genomeSize, args.spadesKmers, maximumReadsLength, args.spadesDefaultKmers, args.spadesMinKmerCovContigs, assembled_se_reads)
					runs['SPAdes'] = [run_successfully, pass_qc, time_taken, failing]

					if run_successfully:
						# Run Pilon
						contigs = contigs_spades

						if not args.skipPilon:
							run_successfully, _, time_taken, failing, assembly_polished, bam_file, pilon_folder = pilon.runPilon(jar_path_pilon, contigs_spades, fastq_files, threads, outdir, jarMaxMemory)
							runs['Pilon'] = [run_successfully, None, time_taken, failing]

							if run_successfully:
								contigs = assembly_polished

							# Run Assembly Mapping check
							if bam_file is not None:
								if not args.skipAssemblyMapping:
									run_successfully, pass_qc, time_taken, failing, assembly_filtered = assembly_mapping.runAssemblyMapping(bam_file, contigs_spades, threads, outdir, args.assemblyMinCoverageContigs, assembly_polished, genomeSize)
									runs['Assembly_Mapping'] = [run_successfully, pass_qc, time_taken, failing]

									if run_successfully:
										contigs = assembly_filtered
								else:
									print '--skipAssemblyMapping set. Skipping Assembly Mapping check'
									runs['Assembly_Mapping'] = skipped
							else:
								print 'Pilon did not produce the bam file! Assembly Mapping check'
								runs['Assembly_Mapping'] = skipped

							if not args.pilonKeepFiles:
								utils.removeDirectory(pilon_folder)

						else:
							print '--skipPilon set. Skipping Pilon correction and Assembly Mapping check'
							runs['Pilon'] = skipped
							runs['Assembly_Mapping'] = skipped

						print '\n' + 'Final assembly: ' + contigs
						with open(os.path.join(outdir, 'final_assembly.txt'), 'wt') as writer:
							writer.write(contigs + '\n')

						# Run MLST
						if not args.skipMLST:
							runs['MLST'] = mlst.runMlst(contigs, scheme, outdir)
						else:
							print '--skipMLST set. Skipping MLST analysis'
							runs['MLST'] = skipped
					else:
						print 'SPAdes did not run successfully! Skipping Pilon correction, Assembly Mapping check and MLST analysis'
						runs['Pilon'] = skipped
						runs['Assembly_Mapping'] = skipped
						runs['MLST'] = skipped

				else:
					print '--skipSPAdes set. Skipping SPAdes, Pilon correction, Assembly Mapping check and MLST analysis'
					runs['SPAdes'] = skipped
					runs['Pilon'] = skipped
					runs['Assembly_Mapping'] = skipped
					runs['MLST'] = skipped
	else:
		print 'Moving to the next sample'
		for step in ('first_Coverage', 'trueCoverage_ReMatCh', 'first_FastQC', 'Trimmomatic', 'second_Coverage', 'second_FastQC', 'Pear', 'SPAdes', 'Pilon', 'Assembly_Mapping', 'MLST'):
			if step == 'Trimmomatic':
				runs[step] = not_run + ['NA']
			else:
				runs[step] = not_run

	# Remove Trimmomatic directory with cleaned reads
	if not args.trimKeepFiles and 'trimmomatic_folder' in locals():
		utils.removeDirectory(trimmomatic_folder)
	# Remove Pear directory
	if not args.pearKeepFiles and 'pear_folder' in locals():
		utils.removeDirectory(pear_folder)

	# Check run
	run_successfully = all(runs[step][0] or runs[step][0] is None for step in runs)

	pass_fastqIntegrity = runs['FastQ_Integrity'][0]
	pass_cov = (runs['second_Coverage'][1] or (runs['second_Coverage'][1] is None and runs['first_Coverage'][1])) is not False
	pass_trueCov = runs['trueCoverage_ReMatCh'][1] is not False
	pass_fastqc = (runs['second_FastQC'][1] or (runs['second_FastQC'][1] is None and runs['first_FastQC'][1])) is not False
	pass_trimmomatic = runs['Trimmomatic'][1] is not False
	pass_spades = runs['Pear'][1] is not False
	pass_spades = runs['SPAdes'][1] is not False
	pass_assemblyMapping = runs['Assembly_Mapping'][1] is not False
	pass_mlst = runs['MLST'][1] is not False
	pass_qc = all([pass_fastqIntegrity, pass_cov, pass_trueCov, pass_fastqc, pass_trimmomatic, pass_spades, pass_assemblyMapping, pass_mlst])

	return run_successfully, pass_qc, runs


if __name__ == "__main__":
	main()
