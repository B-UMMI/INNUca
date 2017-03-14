import utils
import os
from functools import partial
import multiprocessing
import spades


def sequenceHeaders(fastaFile):
	headers = []

	with open(fastaFile, 'rtU') as lines:
		for line in lines:
			line = line.splitlines()[0]
			if len(line) > 0:
				if line[0] == '>':
					headers.append(line[1:])

	if len(headers) == 0:
		headers = None

	return headers


# Get genome coverage data
def compute_genome_coverage_data(alignment_file, sequence_to_analyse, outdir, counter):
	genome_coverage_data_file = os.path.join(outdir, 'samtools_depth.sequence_' + str(counter) + '.tab')
	command = ['samtools', 'depth', '-a', '-r', sequence_to_analyse, alignment_file, '>', genome_coverage_data_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)
	return run_successfully, genome_coverage_data_file


def calculate_genome_coverage(genome_coverage_data_file):
	sequence = None
	position = 0
	coverage = 0
	problems_found = False

	if os.path.getsize(genome_coverage_data_file) > 0:
		with open(genome_coverage_data_file, 'rtU') as reader:
			for line in reader:
				line = line.splitlines()[0]
				if len(line) > 0:
					line = line.split('\t')
					if position == 0 and sequence is None:
						sequence = line[0]
					else:
						if sequence != line[0]:
							print 'WARNING: different sequences found in the same file'
							problems_found = True
							break
					if int(line[1]) == position + 1:
						coverage += int(line[2])
					else:
						print 'WARNING: missing some positions!'
						problems_found = True
						break
					position += 1

	return problems_found, position, coverage


@utils.trace_unhandled_exceptions
def get_sequence_coverage(alignment_file, sequence_to_analyse, outdir, counter):
	problems_found = False
	position = 0
	coverage = 0

	run_successfully, genome_coverage_data_file = compute_genome_coverage_data(alignment_file, sequence_to_analyse, outdir, counter)
	if run_successfully:
		problems_found, position, coverage = calculate_genome_coverage(genome_coverage_data_file)
		if not problems_found and position == 0:
			# Assuming SPAdes headers renamed (with sample name at the beginning)
			position = int(sequence_to_analyse.rsplit('_', 6)[4])
	else:
		problems_found = True

	try:
		os.remove(genome_coverage_data_file)
	except Exception as e:
		print e

	utils.saveVariableToPickle([sequence_to_analyse, run_successfully, problems_found, position, coverage], outdir, str('coverage.sequence_' + str(counter)))


def sample_coverage(referenceFile, alignment_file, outdir, threads):
	coverage_outdir = os.path.join(outdir, 'samtools_depth', '')
	utils.removeDirectory(coverage_outdir)
	os.makedirs(coverage_outdir)

	sequences = sequenceHeaders(referenceFile)

	pool = multiprocessing.Pool(processes=threads)
	counter = 0
	for sequence in sequences:
		pool.apply_async(get_sequence_coverage, args=(alignment_file, sequence, coverage_outdir, counter,))
		counter += 1
	pool.close()
	pool.join()

	sample_coverage_no_problems = True
	mean_coverage_data = {}
	files = [f for f in os.listdir(coverage_outdir) if not f.startswith('.') and os.path.isfile(os.path.join(coverage_outdir, f))]
	for file_found in files:
		if file_found.startswith('coverage.sequence_') and file_found.endswith('.pkl'):
			file_path = os.path.join(coverage_outdir, file_found)

			if sample_coverage_no_problems:
				sequence_to_analyse, run_successfully, problems_found, position, coverage = utils.extractVariableFromPickle(file_path)
				if run_successfully and not problems_found:
					mean_coverage_data[sequence_to_analyse] = {'position': position, 'coverage': coverage, 'mean_coverage': round((float(coverage) / float(position)), 2)}
				else:
					print 'WARNING: it was not possible to compute coverage information for sequence ' + sequence_to_analyse
					sample_coverage_no_problems = False

			os.remove(file_path)

	return sample_coverage_no_problems, mean_coverage_data


def save_assembly_coverage_report(mean_coverage_data, outdir, minCoverageAssembly):
	pass_qc = False
	failing_reason = None

	report_file = os.path.join(outdir, 'assembly_mapping_report.coverage.txt')

	position_initial = 0
	coverage_initial = 0

	# Gather data to determine assembly mean coverage
	for sequence in mean_coverage_data:
		position_initial += mean_coverage_data[sequence]['position']
		coverage_initial += mean_coverage_data[sequence]['coverage']

	# Specify minCoverageAssembly
	if minCoverageAssembly is None:
		if float(coverage_initial) / float(position_initial) / float(3) >= 10:
			minCoverageAssembly = float(coverage_initial) / float(position_initial) / float(3)
			print 'The --assemblyMinCoverageContigs used to filtered low covered contigs was set to one third of the assembly mean coverage (' + str(round((float(coverage_initial) / float(position_initial)), 2)) + 'x): ' + str(round(minCoverageAssembly, 2)) + 'x'
		else:
			minCoverageAssembly = 10
			print 'The --assemblyMinCoverageContigs used to filtered low covered contigs was set to 10 because the assembly mean coverage (' + str(round((float(coverage_initial) / float(position_initial)), 2)) + 'x) was lower than 30x'

	sequences_2_keep = []
	position_filtered = 0
	coverage_filtered = 0
	with open(report_file, 'wt') as writer:
		writer.write('#by_contigs' + '\n')
		for sequence in sorted(mean_coverage_data):
			if mean_coverage_data[sequence]['mean_coverage'] >= minCoverageAssembly:
				sequences_2_keep.append(sequence)
				position_filtered += mean_coverage_data[sequence]['position']
				coverage_filtered += mean_coverage_data[sequence]['coverage']

			writer.write(str('>' + sequence) + '\n' + str(mean_coverage_data[sequence]['mean_coverage']) + '\n')
			writer.flush()
		writer.write('#general' + '\n')
		writer.write('>initial' + '\n' + str(round((float(coverage_initial) / float(position_initial)), 2)) + '\n')
		writer.write('>filtered' + '\n' + str(round((float(coverage_filtered) / float(position_filtered)), 2)) + '\n')

	if round((float(coverage_filtered) / float(position_filtered)), 2) >= 30:
		pass_qc = True
		print 'Assembly coverage: ' + str(round((float(coverage_filtered) / float(position_filtered)), 2)) + 'x'
	else:
		failing_reason = 'Assembly coverage: ' + str(round((float(coverage_filtered) / float(position_filtered)), 2)) + 'x (lower than 30x)'

	return pass_qc, failing_reason, sequences_2_keep


def get_sequence_information(sequence_file):
	sequence_dict = {}

	with open(sequence_file, 'rtU') as original_sequences:
		blank_line_found = False
		sequence_counter = 0
		for line in original_sequences:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not blank_line_found:
					if line.startswith('>'):
						sequence_counter += 1
						sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': [], 'length': 0, 'discard': True}
					else:
						sequence_dict[sequence_counter]['sequence'].append(line)
						sequence_dict[sequence_counter]['length'] += len(line)
				else:
					sequence_dict = None
			else:
				blank_line_found = True

	return sequence_dict


def determine_sequences_to_filter(sequence_dict, list_sequences, pilon_run_successfuly):
	for i in sequence_dict:
		sequence_dict[i]['discard'] = True

	sequence_report_general = {'filtered': {'contigs': 0, 'bp': 0}}

	for i in sequence_dict:
		header = sequence_dict[i]['header']
		if pilon_run_successfuly:
			header = sequence_dict[i]['header'].rsplit('_', 1)[0]

		if header in list_sequences:
			sequence_dict[i]['discard'] = False
			sequence_report_general['filtered']['contigs'] += 1
			sequence_report_general['filtered']['bp'] += sequence_dict[i]['length']

	return sequence_dict, sequence_report_general


def write_filtered_sequences_and_stats(sequence_dict, sequence_report_general, filtered_sequence_file, saveExcludedContigs):
	if saveExcludedContigs:
		path_excluded_contigs = os.path.splitext(filtered_sequence_file)[0] + '.excluded_contigs.fasta'
		excluded_contigs = open(path_excluded_contigs, 'wt')

	found_excluded_contigs = False
	with open(os.path.join(os.path.dirname(filtered_sequence_file), str('assembly_mapping_report.sequences_filtered.' + os.path.splitext(os.path.basename(filtered_sequence_file))[0]) + '.tab'), 'wt') as report_filtered:
		with open(filtered_sequence_file, 'wt') as contigs_filtered:
			fields = ['header', 'length']
			report_filtered.write('\n'.join(['#general', '>contigs', str(sequence_report_general['filtered']['contigs']), '>bp', str(sequence_report_general['filtered']['bp'])]) + '\n')
			report_filtered.write('#' + '\t'.join(fields) + '\n')

			for i in range(1, len(sequence_dict) + 1):
				if not sequence_dict[i]['discard']:
					report_filtered.write('\t'.join([str(sequence_dict[i][f]) for f in fields]) + '\n')
					contigs_filtered.write('>' + sequence_dict[i]['header'] + '\n' + '\n'.join(sequence_dict[i]['sequence']) + '\n')
				else:
					if saveExcludedContigs:
						found_excluded_contigs = True
						excluded_contigs.write('>' + sequence_dict[i]['header'] + '\n' + '\n'.join(sequence_dict[i]['sequence']) + '\n')

	if saveExcludedContigs:
		excluded_contigs.flush()
		excluded_contigs.close()
		if not found_excluded_contigs:
			os.remove(path_excluded_contigs)


def getting_mapping_statistics(alignment_file):
	command = ['samtools', 'flagstat', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)

	dict_mapping_statistics = {}
	if run_successfully:
		stdout = stdout.splitlines()
		for line in stdout:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split(' ', 3)
				field = line[3].split('(', 1)
				if len(field) == 0:
					field = field[0].replace(' ', '_')
				else:
					field = field[0].rsplit(' ', 1)[0].replace(' ', '_')
				dict_mapping_statistics[field] = {'qc_passed': int(line[0]), 'qc_failed': int(line[2])}
	return run_successfully, dict_mapping_statistics


def save_mapping_statistics(dict_mapping_statistics, outdir):
	pass_qc = False
	failing_reason = None

	report_file = os.path.join(outdir, 'assembly_mapping_report.mapping.txt')

	total_reads = 0
	total_mapped_reads = 0
	with open(report_file, 'wt') as writer:
		for field in sorted(dict_mapping_statistics):
			writer.write(str('#' + field) + '\n' + str('>' + 'qc_passed') + '\n' + str(dict_mapping_statistics[field]['qc_passed']) + '\n' + str('>' + 'qc_failed') + '\n' + str(dict_mapping_statistics[field]['qc_failed']) + '\n')

			if field == 'in_total':
				total_reads = dict_mapping_statistics[field]['qc_passed'] + dict_mapping_statistics[field]['qc_failed']
			elif field == 'mapped':
				total_mapped_reads = dict_mapping_statistics[field]['qc_passed'] + dict_mapping_statistics[field]['qc_failed']

	if total_mapped_reads > 0 and total_reads > 0:
		if round((float(total_mapped_reads) / float(total_reads)), 2) >= 0.66:
			pass_qc = True
			print 'Mapped reads: ' + str(round((float(total_mapped_reads) / float(total_reads)), 2) * 100) + '%'
		else:
			failing_reason = 'Mapped reads: ' + str(round((float(total_mapped_reads) / float(total_reads)), 2) * 100) + '% (lower than 66%)'
	else:
		failing_reason = 'No reads were mapped'

	return pass_qc, failing_reason


# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


# Mapping with Bowtie2
def mappingBowtie2(fastq_files, referenceFile, threads, outdir):
	sam_file = os.path.join(outdir, str('alignment.sam'))

	# Index reference file
	run_successfully = indexSequenceBowtie2(referenceFile, threads)

	if run_successfully:
		command = ['bowtie2', '-q', '--very-sensitive-local', '--threads', str(threads), '-x', referenceFile, '', '-S', sam_file]
		if len(fastq_files) == 1:
			command[8] = '-U ' + fastq_files
		else:
			command[8] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

	if not run_successfully:
		sam_file = None

	return run_successfully, sam_file


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
	outFormat_string = os.path.splitext(output_file)[1][1:].lower()
	command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
	if sortByName_True:
		command[6] = '-n'
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	if not run_successfully:
		output_file = None
	return run_successfully, output_file


# Index alignment file
def indexAlignment(alignment_file, print_command_True):
	command = ['samtools', 'index', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, print_command_True)
	return run_successfully


def get_bam_subset(alignment_file, sequences_2_keep, threads):
	bam_subset = os.path.splitext(alignment_file)[0] + '.subset.bam'

	command = ['samtools', 'view', '-buh', '-F', '4', '-o', bam_subset, '-@', str(threads), alignment_file, ' '.join(sequences_2_keep)]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)

	if not run_successfully:
		bam_subset = None

	return run_successfully, bam_subset


assemblyMapping_timer = partial(utils.timer, name='Assembly mapping check')


@assemblyMapping_timer
def runAssemblyMapping(fastq_files, reference_file, threads, outdir, minCoverageAssembly, estimatedGenomeSizeMb, saveExcludedContigs):
	pass_qc = False
	pass_qc_coverage = False
	pass_qc_mapping = False
	pass_qc_sequences = False

	failing = {}

	assemblyMapping_folder = os.path.join(outdir, 'assemblyMapping', '')
	utils.removeDirectory(assemblyMapping_folder)
	os.mkdir(assemblyMapping_folder)

	assembly_filtered = None

	# Create a symbolic link to the assembly
	assembly_link = os.path.join(assemblyMapping_folder, os.path.basename(reference_file))
	os.symlink(reference_file, assembly_link)

	bam_file = None
	# Index assembly using Bowtie2
	run_successfully = indexSequenceBowtie2(assembly_link, threads)

	sample_mapping_statistics_no_problems = False
	if run_successfully:
		run_successfully, sam_file = mappingBowtie2(fastq_files, assembly_link, threads, assemblyMapping_folder)

		if run_successfully:
			bam_file = os.path.splitext(sam_file)[0] + '.bam'
			run_successfully, bam_file = sortAlignment(sam_file, bam_file, False, threads)

			if run_successfully:
				os.remove(sam_file)
				run_successfully = indexAlignment(bam_file, True)

				if run_successfully:
					sequences_2_keep = []
					# Get assembly coverage
					sample_coverage_no_problems, mean_coverage_data = sample_coverage(reference_file, bam_file, assemblyMapping_folder, threads)
					if sample_coverage_no_problems:
						pass_qc_coverage, failing_reason, sequences_2_keep = save_assembly_coverage_report(mean_coverage_data, outdir, minCoverageAssembly)
						if not pass_qc_coverage:
							failing['Coverage'] = [failing_reason]

						assembly_filtered = os.path.splitext(reference_file)[0] + '.mappingCov.fasta'

						sequence_dict, ignore = utils.get_sequence_information(reference_file, 0)
						sequence_dict, sequence_report_general = determine_sequences_to_filter(sequence_dict, sequences_2_keep, False)
						failing_sequences_filtered, minimumBP = spades.qc_assembly(sequence_report_general, estimatedGenomeSizeMb)
						if failing_sequences_filtered['sample'] is not False:
							failing['Sequences_filtered'] = [failing_sequences_filtered['sample']]
							if not minimumBP:
								assembly_filtered = reference_file
							else:
								write_filtered_sequences_and_stats(sequence_dict, sequence_report_general, assembly_filtered, saveExcludedContigs)
						else:
							write_filtered_sequences_and_stats(sequence_dict, sequence_report_general, assembly_filtered, saveExcludedContigs)
							pass_qc_sequences = True
					else:
						failing['Coverage'] = ['Did not run']

					# Save mapping statistics
					sample_mapping_statistics_no_problems, dict_mapping_statistics = getting_mapping_statistics(bam_file)
					if sample_mapping_statistics_no_problems:
						pass_qc_mapping, failing_reason = save_mapping_statistics(dict_mapping_statistics, outdir)

						if not pass_qc_mapping:
							failing['Mapping'] = [failing_reason]
					else:
						failing['Mapping'] = ['Did not run']

					if assembly_filtered is not None and assembly_filtered != reference_file and len(sequences_2_keep) > 0:
						print 'Producing bam subset for sequences to keep'
						run_successfully, bam_subset = get_bam_subset(bam_file, sequences_2_keep, threads)
						if run_successfully:
							os.remove(bam_file)
							os.remove(bam_file + '.bai')
							bam_file = bam_subset
							run_successfully = indexAlignment(bam_file, False)

	if not run_successfully:
		failing['Coverage'] = ['Did not run']

	if len(failing) == 0:
		failing = {'sample': False}
	else:
		print 'Failing:', failing

	run_successfully = sample_coverage_no_problems and sample_mapping_statistics_no_problems
	pass_qc = all([pass_qc_coverage, pass_qc_mapping, pass_qc_sequences])

	if not pass_qc:
		print 'Sample FAILS Assembly Mapping check with: ' + str(failing)

	return run_successfully, pass_qc, failing, assembly_filtered, bam_file, assemblyMapping_folder
