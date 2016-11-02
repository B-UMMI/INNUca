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

	position = 0
	coverage = 0
	sequences_2_keep = []
	with open(report_file, 'wt') as writer:
		writer.write('#by_contigs' + '\n')
		for sequence in sorted(mean_coverage_data):
			position += mean_coverage_data[sequence]['position']
			coverage += mean_coverage_data[sequence]['coverage']

			if mean_coverage_data[sequence]['mean_coverage'] >= minCoverageAssembly:
				sequences_2_keep.append(sequence)

			writer.write(str('>' + sequence) + '\n' + str(mean_coverage_data[sequence]['mean_coverage']) + '\n')
			writer.flush()
		writer.write('#general' + '\n' + str(round((float(coverage) / float(position)), 2)) + '\n')

	if round((float(coverage) / float(position)), 2) >= 30:
		pass_qc = True
		print 'Assembly coverage: ' + str(round((float(coverage) / float(position)), 2)) + 'x'
	else:
		failing_reason = 'Assembly coverage: ' + str(round((float(coverage) / float(position)), 2)) + 'x (lower than 30x)'

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


def write_filtered_sequences_and_stats(sequence_dict, sequence_report_general, filtered_sequence_file):
	with open(os.path.join(os.path.dirname(filtered_sequence_file), str('assembly_mapping_report.sequences_filtered.' + os.path.splitext(os.path.basename(filtered_sequence_file))[0]) + '.tab'), 'wt') as report_filtered:
		with open(filtered_sequence_file, 'wt') as contigs_filtered:
			fields = ['header', 'length']
			report_filtered.write('\n'.join(['#general', '>contigs', str(sequence_report_general['filtered']['contigs']), '>bp', str(sequence_report_general['filtered']['bp'])]) + '\n')
			report_filtered.write('#' + '\t'.join(fields) + '\n')

			for i in range(1, len(sequence_dict) + 1):
				if not sequence_dict[i]['discard']:
					report_filtered.write('\t'.join([str(sequence_dict[i][f]) for f in fields]) + '\n')
					contigs_filtered.write('>' + sequence_dict[i]['header'] + '\n' + '\n'.join(sequence_dict[i]['sequence']) + '\n')


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


assemblyMapping_timer = partial(utils.timer, name='Assembly mapping check')


@assemblyMapping_timer
def runAssemblyMapping(alignment_file, reference_file, threads, outdir, minCoverageAssembly, assembly_pilon, keep_pilon_assembly, estimatedGenomeSizeMb):
	pass_qc = False
	pass_qc_coverage = False
	pass_qc_mapping = False
	pass_qc_sequences = False

	failing = {}

	assemblyMapping_folder = os.path.join(outdir, 'assemblyMapping', '')
	utils.removeDirectory(assemblyMapping_folder)
	os.mkdir(assemblyMapping_folder)

	assembly_filtered = None

	pilon_run_successfuly = True if assembly_pilon is not None else False

	# Get assembly coverage
	sample_coverage_no_problems, mean_coverage_data = sample_coverage(reference_file, alignment_file, assemblyMapping_folder, threads)
	if sample_coverage_no_problems:
		pass_qc_coverage, failing_reason, sequences_2_keep = save_assembly_coverage_report(mean_coverage_data, outdir, minCoverageAssembly)
		if not pass_qc_coverage:
			failing['Coverage'] = [failing_reason]

		assembly = reference_file if assembly_pilon is None else assembly_pilon
		assembly_filtered = os.path.splitext(assembly)[0] + '.mappingCov.fasta'

		sequence_dict = get_sequence_information(assembly)
		sequence_dict, sequence_report_general = determine_sequences_to_filter(sequence_dict, sequences_2_keep, pilon_run_successfuly)
		failing_sequences_filtered = spades.qc_assembly(sequence_report_general, estimatedGenomeSizeMb)
		if failing_sequences_filtered['sample'] is False:
			write_filtered_sequences_and_stats(sequence_dict, sequence_report_general, assembly_filtered)
			pass_qc_sequences = True
		else:
			failing['Sequences_filtered'] = [failing_sequences_filtered['sample']]
			assembly_filtered = assembly
	else:
		failing['Coverage'] = ['Did not run']

	if pilon_run_successfuly and not keep_pilon_assembly:
		os.remove(assembly_pilon)

	# Save mapping statistics
	sample_mapping_statistics_no_problems, dict_mapping_statistics = getting_mapping_statistics(alignment_file)
	if sample_mapping_statistics_no_problems:
		pass_qc_mapping, failing_reason = save_mapping_statistics(dict_mapping_statistics, outdir)

		if not pass_qc_mapping:
			failing['Mapping'] = [failing_reason]
	else:
		failing['Mapping'] = ['Did not run']

	run_successfully = sample_coverage_no_problems and sample_mapping_statistics_no_problems
	pass_qc = all([pass_qc_coverage, pass_qc_mapping, pass_qc_sequences])

	if not pass_qc:
		print 'Sample FAILS Assembly Mapping check with: ' + str(failing)

	utils.removeDirectory(assemblyMapping_folder)

	return run_successfully, pass_qc, failing, assembly_filtered
