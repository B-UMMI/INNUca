import pilon
import assembly_mapping
import utils
import os
import multiprocessing
import sys
import functools


def check_existing_default_config(species, script_path):
	species = species.lower().split(' ')
	trueCoverage_config_folder = os.path.join(os.path.dirname(script_path), 'modules', 'trueCoverage_rematch', '')
	config = None
	reference = None
	files = [f for f in os.listdir(trueCoverage_config_folder) if not f.startswith('.') and os.path.isfile(os.path.join(trueCoverage_config_folder, f))]
	for file_found in files:
		file_path = os.path.join(trueCoverage_config_folder, file_found)
		if file_found == '_'.join(species) + '.config':
			config = file_path
		elif file_found == '_'.join(species) + '.fasta':
			reference = file_path
	return config, reference


def parse_config(config_file):
	config = {'reference_file': None, 'length_extra_seq': None, 'maximum_number_absent_genes': None, 'maximum_number_genes_multiple_alleles': None, 'minimum_read_coverage': None, 'minimum_depth_presence': None, 'minimum_depth_call': None, 'minimum_depth_frequency_dominant_allele': None, 'minimum_gene_coverage': None}

	with open(config_file, 'rtU') as reader:
		field = None
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split(' ')[0]
				if line.startswith('#'):
					line = line[1:].split(' ')[0]
					field = line
				else:
					if field is not None:
						if field in ['length_extra_seq', 'maximum_number_absent_genes', 'maximum_number_genes_multiple_alleles', 'minimum_read_coverage', 'minimum_depth_presence', 'minimum_depth_call', 'minimum_gene_coverage']:
							line = int(line)
							if field == 'minimum_gene_coverage':
								if line < 0 or line > 100:
									sys.exit('minimum_gene_coverage in trueCoverage_rematch config file must be an integer between 0 and 100')
						elif field == 'minimum_depth_frequency_dominant_allele':
							line = float(line)
							if line < 0 or line > 1:
								sys.exit('minimum_depth_frequency_dominant_allele in trueCoverage_rematch config file must be a double between 0 and 1')
						config[field] = line
						field = None

	for field in config:
		if config[field] is None:
			sys.exit(field + ' in trueCoverage_rematch config file is missing')

	return config


def index_fasta_samtools(fasta):
	command = ['samtools', 'faidx', fasta]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


def mapping_reads(fastq_files, reference_file, threads, outdir):
	# Create a symbolic link to the reference_file
	reference_link = os.path.join(outdir, os.path.basename(reference_file))
	os.symlink(reference_file, reference_link)

	# Index assembly using Bowtie2
	run_successfully = pilon.indexSequenceBowtie2(reference_link, threads)

	bam_file = None
	if run_successfully:
		# Mapping reads using Bowtie2
		run_successfully, sam_file = pilon.mappingBowtie2(fastq_files, reference_link, threads, outdir)

		if run_successfully:
			# Convert sam to bam and sort bam
			run_successfully, bam_file = pilon.sortAlignment(sam_file, str(os.path.splitext(sam_file)[0] + '.bam'), False, threads)

			if run_successfully:
				os.remove(sam_file)
				# Index bam
				run_successfully = pilon.indexAlignment(bam_file)

	return run_successfully, bam_file, reference_link


def create_vcf(bam_file, sequence_to_analyse, outdir, counter, reference_file):
	gene_vcf = os.path.join(outdir, 'samtools_mpileup.sequence_' + str(counter) + '.vcf')

	command = ['samtools', 'mpileup', '--count-orphans', '--no-BAQ', '--min-BQ', '0', '--min-MQ', '0', '--fasta-ref', reference_file, '--region', sequence_to_analyse, '--output', gene_vcf, '--VCF', '--uncompressed', '--output-tags', 'INFO/AD,AD,DP', bam_file]

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
	if not run_successfully:
		gene_vcf = None
	return run_successfully, gene_vcf


# Read vcf file
class Vcf():
	def __init__(self, vcfFile):
		self.vcf = open(vcfFile, 'rtU')
		self.line_read = self.vcf.readline()
		while self.line_read.startswith('#'):
			self.line_read = self.vcf.readline()
		self.line = self.line_read

	def readline(self):
		self.line_stored = self.line
		self.line = self.vcf.readline()
		return self.line_stored

	def close(self):
		self.vcf.close()


def get_variants(gene_vcf):
	variants = {}

	vfc_file = Vcf(gene_vcf)
	line = vfc_file.readline()
	while len(line) > 0:
		fields = line.splitlines()[0].split('\t')
		if len(fields) > 0:
			fields[1] = int(fields[1])

			info_field = {}
			for i in fields[7].split(';'):
				i = i.split('=')
				if len(i) > 1:
					info_field[i[0]] = i[1]
				else:
					info_field[i[0]] = None

			format_field = {}
			format_field_name = fields[8].split(':')
			format_data = fields[9].split(':')

			for i in range(0, len(format_data)):
				format_field[format_field_name[i]] = format_data[i].split(',')

			fields_to_store = {'REF': fields[3], 'ALT': fields[4].split(','), 'info': info_field, 'format': format_field}
			if fields[1] in variants:
				variants[fields[1]][len(variants[fields[1]])] = fields_to_store
			else:
				variants[fields[1]] = {0: fields_to_store}

		line = vfc_file.readline()
	vfc_file.close()

	return variants


def indel_entry(variant_position):
	entry_with_indel = None
	for i in variant_position:
		keys = variant_position[i]['info'].keys()
		if 'INDEL' in keys:
			entry_with_indel = i

	return entry_with_indel


def determine_variant(variant_position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	alt = None
	low_coverage = False
	multiple_alleles = False

	if int(variant_position['format']['DP'][0]) >= minimum_depth_presence:
		if int(variant_position['format']['DP'][0]) < minimum_depth_call:
			alt = 'N'
			low_coverage = True
		else:
			index_dominant_allele = variant_position['format']['AD'].index(str(max(map(int, variant_position['format']['AD']))))

			if int(variant_position['format']['AD'][index_dominant_allele]) < minimum_depth_call:
				alt = 'N'
				low_coverage = True
			else:
				if float(variant_position['format']['AD'][index_dominant_allele]) / float(variant_position['format']['DP'][0]) < minimum_depth_frequency_dominant_allele:
					multiple_alleles = True
					alt = 'N'
				else:
					if index_dominant_allele == 0:
						alt = '.'
					else:
						alt = variant_position['ALT'][index_dominant_allele - 1]
	else:
		low_coverage = True

	return variant_position['REF'], alt, low_coverage, multiple_alleles


def snp_indel(variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	ref = alt = low_coverage = None
	multiple_alleles = False

	entry_with_indel = indel_entry(variants[position])
	if entry_with_indel is None:
		ref, alt, low_coverage, multiple_alleles = determine_variant(variants[position][0], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)
	else:
		ref, alt, low_coverage, multiple_alleles = determine_variant(variants[position][entry_with_indel], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)
		if alt is not None and (low_coverage or multiple_alleles):
			alt = 'N' * len(ref)

	return ref, alt, low_coverage, multiple_alleles


def find_multiple_alleles(variants, sequence_length, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, length_extra_seq):
	multiple_alleles_found = 0

	counter = 1
	while counter <= sequence_length:
		if counter in variants:
			ref, alt, low_coverage, multiple_alleles = snp_indel(variants, counter, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)
			if counter <= length_extra_seq and counter + len(ref) > length_extra_seq:
				if multiple_alleles:
					multiple_alleles_found += 1
			elif counter > length_extra_seq and counter <= sequence_length - length_extra_seq:
				if multiple_alleles:
					multiple_alleles_found += 1
			counter += len(ref)
		else:
			counter += 1

	return multiple_alleles_found


def get_coverage(gene_coverage):
	coverage = {}

	with open(gene_coverage, 'rtU') as reader:
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split('\t')
				coverage[int(line[1])] = int(line[2])

	return coverage


def get_coverage_report(coverage, sequence_length, minimum_depth_presence, minimum_depth_call, length_extra_seq):
	if len(coverage) == 0:
		return 100.0, 100.0, 0.0

	count_absent = 0
	count_lowCoverage = 0
	sum_coverage = 0

	counter = 1
	while counter <= sequence_length:
		if counter > length_extra_seq and counter <= sequence_length - length_extra_seq:
			if coverage[counter] < minimum_depth_presence:
				count_absent += 1
				count_lowCoverage += 1
			else:
				if coverage[counter] < minimum_depth_call:
					count_lowCoverage += 1
				sum_coverage += coverage[counter]
		counter += 1

	mean_coverage = 0
	if sequence_length - 2 * length_extra_seq - count_absent > 0:
		mean_coverage = float(sum_coverage) / float(sequence_length - 2 * length_extra_seq - count_absent)

	return float(count_absent) / float(sequence_length - 2 * length_extra_seq) * 100, float(count_lowCoverage) / float(sequence_length - 2 * length_extra_seq) * 100, mean_coverage


@utils.trace_unhandled_exceptions
def analyse_sequence_data(bam_file, sequence_information, outdir, counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	multiple_alleles_found = percentage_absent = percentage_lowCoverage = meanCoverage = None

	# Create vcf file (for multiple alleles check)
	run_successfully, gene_vcf = create_vcf(bam_file, sequence_information['header'], outdir, counter, reference_file)

	if run_successfully:
		# Create coverage tab file
		run_successfully, gene_coverage = assembly_mapping.compute_genome_coverage_data(bam_file, sequence_information['header'], outdir, counter)

		if run_successfully:
			variants = get_variants(gene_vcf)
			coverage = get_coverage(gene_coverage)

			multiple_alleles_found = find_multiple_alleles(variants, sequence_information['length'], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, length_extra_seq)

			percentage_absent, percentage_lowCoverage, meanCoverage = get_coverage_report(coverage, sequence_information['length'], minimum_depth_presence, minimum_depth_call, length_extra_seq)

	utils.saveVariableToPickle([run_successfully, counter, multiple_alleles_found, percentage_absent, percentage_lowCoverage, meanCoverage], outdir, str('coverage_info.' + str(counter)))


def get_sequence_information(fasta_file):
	sequence_dict = {}

	with open(fasta_file, 'rtU') as reader:
		blank_line_found = False
		sequence_counter = 0
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not blank_line_found:
					if line.startswith('>'):
						sequence_counter += 1
						sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': [], 'length': 0}
					else:
						sequence_dict[sequence_counter]['sequence'].append(line)
						sequence_dict[sequence_counter]['length'] += len(line)
				else:
					sequence_dict = None
			else:
				blank_line_found = True

	return sequence_dict


def sequence_data(reference_file, bam_file, outdir, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	sequence_data_outdir = os.path.join(outdir, 'sequence_data', '')
	utils.removeDirectory(sequence_data_outdir)
	os.makedirs(sequence_data_outdir)

	sequences = get_sequence_information(reference_file)

	pool = multiprocessing.Pool(processes=threads)
	for sequence_counter in sequences:
		sequence_dir = os.path.join(sequence_data_outdir, str(sequence_counter), '')
		utils.removeDirectory(sequence_dir)
		os.makedirs(sequence_dir)
		pool.apply_async(analyse_sequence_data, args=(bam_file, sequences[sequence_counter], sequence_dir, sequence_counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele,))
	pool.close()
	pool.join()

	run_successfully, sample_data = gather_gene_data_together(sequence_data_outdir, sequences)

	return run_successfully, sample_data


def gather_gene_data_together(data_directory, sequences_information):
	run_successfully = True
	counter = 0
	sample_data = {}

	genes_directories = [d for d in os.listdir(data_directory) if not d.startswith('.') and os.path.isdir(os.path.join(data_directory, d, ''))]
	for gene_dir in genes_directories:
		gene_dir_path = os.path.join(data_directory, gene_dir, '')

		files = [f for f in os.listdir(gene_dir_path) if not f.startswith('.') and os.path.isfile(os.path.join(gene_dir_path, f))]
		for file_found in files:
			if file_found.startswith('coverage_info.') and file_found.endswith('.pkl'):
				file_path = os.path.join(gene_dir_path, file_found)

				if run_successfully:
					run_successfully, sequence_counter, multiple_alleles_found, percentage_absent, percentage_lowCoverage, meanCoverage = utils.extractVariableFromPickle(file_path)
					sample_data[sequence_counter] = {'header': sequences_information[sequence_counter]['header'], 'gene_coverage': 100 - percentage_absent, 'gene_low_coverage': percentage_lowCoverage, 'gene_number_positions_multiple_alleles': multiple_alleles_found, 'gene_mean_read_coverage': meanCoverage}
					counter += 1

		utils.removeDirectory(gene_dir_path)

	if counter != len(sequences_information):
		run_successfully = False

	return run_successfully, sample_data


trueCoverage_timer = functools.partial(utils.timer, name='True coverage check')


@trueCoverage_timer
def runTrueCoverage(fastq_files, reference_file, threads, outdir, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, minimum_gene_coverage, maximum_number_absent_genes, maximum_number_genes_multiple_alleles, minimum_read_coverage):
	pass_qc = False
	failing = {}

	trueCoverage_folder = os.path.join(outdir, 'trueCoverage', '')
	utils.removeDirectory(trueCoverage_folder)
	os.mkdir(trueCoverage_folder)

	# Map reads
	run_successfully, bam_file, reference_file = mapping_reads(fastq_files, reference_file, threads, trueCoverage_folder)

	if run_successfully:
		# Index reference file
		run_successfully = index_fasta_samtools(reference_file)
		if run_successfully:
			run_successfully, sample_data = sequence_data(reference_file, bam_file, trueCoverage_folder, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

			if run_successfully:
				number_absent_genes = 0
				number_genes_multiple_alleles = 0
				mean_sample_coverage = 0
				with open(os.path.join(outdir, 'trueCoverage_report.txt'), 'wt') as writer:
					writer.write('\t'.join(['#gene', 'percentage_gene_coverage', 'gene_mean_read_coverage', 'percentage_gene_low_coverage', 'number_positions_multiple_alleles']) + '\n')
					for i in range(1, len(sample_data) + 1):
						writer.write('\t'.join([sample_data[i]['header'], str(round(sample_data[i]['gene_coverage'], 2)), str(round(sample_data[i]['gene_mean_read_coverage'], 2)), str(round(sample_data[i]['gene_low_coverage'], 2)), str(sample_data[i]['gene_number_positions_multiple_alleles'])]) + '\n')

						if sample_data[i]['gene_coverage'] < minimum_gene_coverage:
							number_absent_genes += 1
						else:
							mean_sample_coverage += sample_data[i]['gene_mean_read_coverage']
							if sample_data[i]['gene_number_positions_multiple_alleles'] > 0:
								number_genes_multiple_alleles += 1

					mean_sample_coverage = float(mean_sample_coverage) / float(len(sample_data) - number_absent_genes)
					writer.write('\n'.join(['#general', '>number_absent_genes', str(number_absent_genes), '>number_genes_multiple_alleles', str(number_genes_multiple_alleles), '>mean_sample_coverage', str(round(mean_sample_coverage, 2))]) + '\n')

					print '\n'.join([str('number_absent_genes: ' + str(number_absent_genes)), str('number_genes_multiple_alleles: ' + str(number_genes_multiple_alleles)), str('mean_sample_coverage: ' + str(round(mean_sample_coverage, 2)))])

				if number_absent_genes > maximum_number_absent_genes:
					failing['absent_genes'] = 'The number of absent genes (' + str(number_absent_genes) + ') exceeds the maximum allowed (' + str(maximum_number_absent_genes) + ')'
				if number_genes_multiple_alleles > maximum_number_genes_multiple_alleles:
					failing['multiple_alleles'] = 'The number of genes with multiple alleles (' + str(number_genes_multiple_alleles) + ') exceeds the maximum allowed (' + str(maximum_number_genes_multiple_alleles) + ')'
				if mean_sample_coverage < minimum_read_coverage:
					failing['read_coverage'] = 'The mean read coverage for genes present (' + str(mean_sample_coverage) + ') dit not meet the minimum required (' + str(minimum_read_coverage) + ')'
			else:
				failing['sample'] = 'Did not run'
		else:
			failing['sample'] = 'Did not run'
	else:
		failing['sample'] = 'Did not run'

	if len(failing) == 0:
		pass_qc = True
		failing['sample'] = False
	else:
		print failing

	utils.removeDirectory(trueCoverage_folder)

	return run_successfully, pass_qc, failing
