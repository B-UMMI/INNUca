import os.path
import multiprocessing
import utils
import functools
import sys


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


def index_fasta_samtools(fasta, region_None, region_outfile_none, print_comand_True):
	command = ['samtools', 'faidx', fasta, '', '', '']
	shell_true = False
	if region_None is not None:
		command[3] = region_None
	if region_outfile_none is not None:
		command[4] = '>'
		command[5] = region_outfile_none
		shell_true = True
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, shell_true, None, print_comand_True)
	return run_successfully, stdout


# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


# Mapping with Bowtie2
def mappingBowtie2(fastq_files, referenceFile, threads, outdir, conserved_True):
	sam_file = os.path.join(outdir, str('alignment.sam'))

	# Index reference file
	run_successfully = indexSequenceBowtie2(referenceFile, threads)

	if run_successfully:
		command = ['bowtie2', '-q', '', '--threads', str(threads), '-x', referenceFile, '', '--no-unal', '-S', sam_file]

		if len(fastq_files) == 1:
			command[7] = '-U ' + fastq_files
		elif len(fastq_files) == 2:
			command[7] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]
		else:
			return False, None

		if conserved_True:
			command[2] = '--sensitive'
		else:
			command[2] = '--very-sensitive-local'

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
def indexAlignment(alignment_file):
	command = ['samtools', 'index', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


def mapping_reads(fastq_files, reference_file, threads, outdir, conserved_True):
	# Create a symbolic link to the reference_file
	reference_link = os.path.join(outdir, os.path.basename(reference_file))
	os.symlink(reference_file, reference_link)

	bam_file = None
	# Mapping reads using Bowtie2
	run_successfully, sam_file = mappingBowtie2(fastq_files, reference_link, threads, outdir, conserved_True)

	if run_successfully:
		# Convert sam to bam and sort bam
		run_successfully, bam_file = sortAlignment(sam_file, str(os.path.splitext(sam_file)[0] + '.bam'), False, threads)

		if run_successfully:
			os.remove(sam_file)
			# Index bam
			run_successfully = indexAlignment(bam_file)

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
	entry_with_indel = []
	entry_with_snp = None
	for i in variant_position:
		keys = variant_position[i]['info'].keys()
		if 'INDEL' in keys:
			entry_with_indel.append(i)
		else:
			entry_with_snp = i

	return entry_with_indel, entry_with_snp


def get_alt_noMatter(variant_position, indel_true):
	# dp = int(variant_position['info']['DP'])
	dp = sum(map(int, variant_position['format']['AD']))
	index_dominant_allele = None
	if not indel_true:
		ad_idv = max(map(int, variant_position['format']['AD']))

		index_dominant_allele = variant_position['format']['AD'].index(str(ad_idv))

		if index_dominant_allele == 0:
			alt = '.'
		else:
			alt = variant_position['ALT'][index_dominant_allele - 1]
	else:
		if float(variant_position['info']['IDV']) / float(dp) > 0.5:
			alt = variant_position['ALT'][0]
			ad_idv = int(variant_position['info']['IDV'])
		else:
			alt = '.'
			ad_idv = max(map(int, variant_position['format']['AD']))

	return alt, dp, ad_idv, index_dominant_allele


def get_alt_correct(variant_position, alt_noMatter, dp, ad_idv, index_dominant_allele, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	alt = None
	low_coverage = False
	multiple_alleles = False

	if dp >= minimum_depth_presence:
		if dp < minimum_depth_call:
			alt = 'N' * len(variant_position['REF'])
			low_coverage = True
		else:
			if ad_idv < minimum_depth_call:
				alt = 'N' * len(variant_position['REF'])
				low_coverage = True
			else:
				if float(ad_idv) / float(dp) < minimum_depth_frequency_dominant_allele:
					if index_dominant_allele is not None:
						variants_coverage = [int(variant_position['format']['AD'][i]) for i in range(0, len([j for j in variant_position['ALT'] if j != '<*>']) + 1) if i != index_dominant_allele]
						if float(max(variants_coverage)) / float(sum(variants_coverage)) > 0.5:
							multiple_alleles = True
							alt = 'N' * len(variant_position['REF'])
						elif float(max(variants_coverage)) / float(sum(variants_coverage)) == 0.5 and len(variants_coverage) > 2:
							multiple_alleles = True
							alt = 'N' * len(variant_position['REF'])
						else:
							alt = alt_noMatter
					else:
						multiple_alleles = True
						alt = 'N' * len(variant_position['REF'])
				else:
					alt = alt_noMatter
	else:
		low_coverage = True

	return alt, low_coverage, multiple_alleles


def get_alt_alignment(ref, alt):
	if alt is None:
		alt = 'N'
	else:
		if len(ref) != len(alt):
			if len(alt) < len(ref):
				alt += 'N' * (len(ref) - len(alt))
			else:
				if alt[:len(ref)] == ref:
					alt = '.'
				else:
					alt = alt[:len(ref)]

	return alt


def get_indel_more_likely(variant_position, indels_entry):
	indel_coverage = {}
	for i in indels_entry:
		indel_coverage[i] = int(variant_position['info']['IDV'])
	return indel_coverage.index(str(max(indel_coverage.values())))


def determine_variant(variant_position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, indel_true):
	alt_noMatter, dp, ad_idv, index_dominant_allele = get_alt_noMatter(variant_position, indel_true)

	alt_correct, low_coverage, multiple_alleles = get_alt_correct(variant_position, alt_noMatter, dp, ad_idv, index_dominant_allele, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

	alt_alignment = get_alt_alignment(variant_position['REF'], alt_correct)

	return variant_position['REF'], alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment


def snp_indel(variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	entry_with_indel, entry_with_snp = indel_entry(variants[position])

	if len(entry_with_indel) == 0:
		ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
	else:
		ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
		if alt_correct is not None and not low_coverage:
			alt_correct_snp = alt_correct
			alt_noMatter_snp = alt_noMatter
			alt_alignment_snp = alt_alignment
			if len(entry_with_indel) > 1:
				indel_more_likely = get_indel_more_likely(variants[position], entry_with_indel)
				ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_indel[indel_more_likely]], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, True)
			else:
				ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_indel[0]], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, True)
			if alt_correct is not None and not low_coverage:
				if alt_correct == '.':
					ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
				else:
					if alt_correct_snp != '.' and alt_correct[0] != alt_correct_snp:
						alt_correct = alt_correct_snp + alt_correct[1:] if len(alt_correct) > 1 else alt_correct_snp
					if alt_noMatter_snp != '.' and alt_noMatter[0] != alt_noMatter_snp:
						alt_noMatter = alt_noMatter_snp + alt_noMatter[1:] if len(alt_noMatter) > 1 else alt_noMatter_snp
					if alt_alignment_snp != '.' and alt_alignment[0] != alt_alignment_snp:
						alt_alignment = alt_alignment_snp + alt_alignment[1:] if len(alt_alignment) > 1 else alt_alignment_snp
			else:
				ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)

	return ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment


def get_true_variants(variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence):
	variants_correct = {}
	variants_noMatter = {}
	variants_alignment = {}

	absent_positions = {}
	last_absent_position = ''
	multiple_alleles_found = []

	counter = 1
	while counter <= len(sequence):
		if counter in variants:
			ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = snp_indel(variants, counter, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

			if alt_alignment != '.':
				variants_alignment[counter] = {'REF': ref, 'ALT': alt_alignment}

			if alt_correct is None:
				if counter - len(last_absent_position) in absent_positions:
					absent_positions[counter - len(last_absent_position)]['REF'] += ref
				else:
					absent_positions[counter] = {'REF': ref, 'ALT': ''}
				last_absent_position += ref
			else:
				if alt_noMatter != '.':
					variants_noMatter[counter] = {'REF': ref, 'ALT': alt_noMatter}

				if alt_correct != '.':
					if len(alt_correct) < len(ref):
						if len(alt_correct) == 1:
							absent_positions[counter + 1] = {'REF': ref[1:], 'ALT': ''}
						else:
							absent_positions[counter + 1] = {'REF': ref[1:], 'ALT': alt_correct[1:]}

						last_absent_position = ref[1:]
					else:
						variants_correct[counter] = {'REF': ref, 'ALT': alt_correct}
						last_absent_position = ''
				else:
					last_absent_position = ''

			if multiple_alleles:
				multiple_alleles_found.append(counter)

			counter += len(ref)
		else:
			variants_alignment[counter] = {'REF': sequence[counter - 1], 'ALT': 'N'}

			if counter - len(last_absent_position) in absent_positions:
				absent_positions[counter - len(last_absent_position)]['REF'] += sequence[counter - 1]
			else:
				absent_positions[counter] = {'REF': sequence[counter - 1], 'ALT': ''}
			last_absent_position += sequence[counter - 1]
			counter += 1

	for position in absent_positions:
		if position == 1:
			variants_correct[position] = {'REF': absent_positions[position]['REF'], 'ALT': 'N'}
			variants_noMatter[position] = {'REF': absent_positions[position]['REF'], 'ALT': 'N'}
		else:
			if position - 1 not in variants_correct:
				variants_correct[position - 1] = {'REF': sequence[position - 2] + absent_positions[position]['REF'], 'ALT': sequence[position - 2] + absent_positions[position]['ALT']}
			else:
				variants_correct[position - 1] = {'REF': variants_correct[position - 1]['REF'] + absent_positions[position]['REF'][len(variants_correct[position - 1]['REF']) - 1:], 'ALT': variants_correct[position - 1]['ALT'] + absent_positions[position]['ALT'][len(variants_correct[position - 1]['ALT']) - 1 if len(variants_correct[position - 1]['ALT']) > 0 else 0:]}

			if position - 1 not in variants_noMatter:
				variants_noMatter[position - 1] = {'REF': sequence[position - 2] + absent_positions[position]['REF'], 'ALT': sequence[position - 2] + absent_positions[position]['ALT']}
			else:
				variants_noMatter[position - 1] = {'REF': variants_noMatter[position - 1]['REF'] + absent_positions[position]['REF'][len(variants_noMatter[position - 1]['REF']) - 1:], 'ALT': variants_noMatter[position - 1]['ALT'] + absent_positions[position]['ALT'][len(variants_noMatter[position - 1]['ALT']) - 1 if len(variants_noMatter[position - 1]['ALT']) > 0 else 0:]}

	return variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found


def clean_variant_in_extra_seq_left(variant_dict, position, length_extra_seq, multiple_alleles_found, number_multi_alleles):
	if position + len(variant_dict[position]['REF']) - 1 > length_extra_seq:
		if multiple_alleles_found is not None and position in multiple_alleles_found:
			number_multi_alleles += 1

		temp_variant = variant_dict[position]
		del variant_dict[position]
		variant_dict[length_extra_seq + 1] = {}
		variant_dict[length_extra_seq + 1]['REF'] = temp_variant['REF'][length_extra_seq - position + 1:]
		variant_dict[length_extra_seq + 1]['ALT'] = temp_variant['ALT'][length_extra_seq - position + 1:] if len(temp_variant['ALT']) >= length_extra_seq - position + 1 else 'N'
	else:
		del variant_dict[position]

	return variant_dict, number_multi_alleles


def clean_variant_in_extra_seq_rigth(variant_dict, position, sequence_length, length_extra_seq):
	if position + len(variant_dict[position]['REF']) - 1 > sequence_length - length_extra_seq:
		variant_dict[position]['REF'] = variant_dict[position]['REF'][: - (position - (sequence_length - length_extra_seq)) + 1]
		variant_dict[position]['ALT'] = variant_dict[position]['ALT'][: - (position - (sequence_length - length_extra_seq)) + 1] if len(variant_dict[position]['ALT']) >= - (position - (sequence_length - length_extra_seq)) + 1 else variant_dict[position]['ALT']

	return variant_dict


def cleanning_variants_extra_seq(variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found, length_extra_seq, sequence_length):
	number_multi_alleles = 0

	counter = 1
	while counter <= sequence_length:
		if counter <= length_extra_seq:
			if counter in variants_correct:
				variants_correct, number_multi_alleles = clean_variant_in_extra_seq_left(variants_correct, counter, length_extra_seq, multiple_alleles_found, number_multi_alleles)
			if counter in variants_noMatter:
				variants_noMatter, ignore = clean_variant_in_extra_seq_left(variants_noMatter, counter, length_extra_seq, None, None)
			if counter in variants_alignment:
				variants_alignment, ignore = clean_variant_in_extra_seq_left(variants_alignment, counter, length_extra_seq, None, None)
		elif counter > length_extra_seq and counter <= sequence_length - length_extra_seq:
			if counter in variants_correct:
				if counter in multiple_alleles_found:
					number_multi_alleles += 1

				variants_correct = clean_variant_in_extra_seq_rigth(variants_correct, counter, sequence_length, length_extra_seq)
			if counter in variants_noMatter:
				variants_noMatter = clean_variant_in_extra_seq_rigth(variants_noMatter, counter, sequence_length, length_extra_seq)
			if counter in variants_alignment:
				variants_alignment = clean_variant_in_extra_seq_rigth(variants_alignment, counter, sequence_length, length_extra_seq)
		else:
			if counter in variants_correct:
				del variants_correct[counter]
			if counter in variants_noMatter:
				del variants_noMatter[counter]
			if counter in variants_alignment:
				del variants_alignment[counter]

		counter += 1

	return variants_correct, variants_noMatter, variants_alignment, number_multi_alleles


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


# Get genome coverage data
def compute_genome_coverage_data(alignment_file, sequence_to_analyse, outdir, counter):
	genome_coverage_data_file = os.path.join(outdir, 'samtools_depth.sequence_' + str(counter) + '.tab')
	command = ['samtools', 'depth', '-a', '-q', '0', '-r', sequence_to_analyse, alignment_file, '>', genome_coverage_data_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)
	return run_successfully, genome_coverage_data_file


# # NOT_USED
# def write_variants_vcf(variants, outdir, sequence_to_analyse, sufix):
# 	vcf_file = os.path.join(outdir, str(sequence_to_analyse + '.' + sufix + '.vcf'))
# 	with open(vcf_file, 'wt') as writer:
# 		writer.write('##fileformat=VCFv4.2' + '\n')
# 		writer.write('#' + '\t'.join(['SEQUENCE', 'POSITION', 'ID_unused', 'REFERENCE_sequence', 'ALTERNATIVE_sequence', 'QUALITY_unused', 'FILTER_unused', 'INFO_unused', 'FORMAT_unused']) + '\n')
# 		for i in sorted(variants.keys()):
# 			writer.write('\t'.join([sequence_to_analyse, str(i), '.', variants[i]['REF'], variants[i]['ALT'], '.', '.', '.', '.']) + '\n')
#
# 	compressed_vcf_file = vcf_file + '.gz'
# 	command = ['bcftools', 'convert', '-o', compressed_vcf_file, '-O', 'z', vcf_file]
# 	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
# 	if run_successfully:
# 		command = ['bcftools', 'index', compressed_vcf_file]
# 		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
#
# 	if not run_successfully:
# 		compressed_vcf_file = None
#
# 	return run_successfully, compressed_vcf_file
#
#
# # NOT_USED
# def parse_fasta_inMemory(fasta_memory):
# 	fasta_memory = fasta_memory.splitlines()
# 	sequence_dict = {}
# 	for line in fasta_memory:
# 		if len(line) > 0:
# 			if line.startswith('>'):
# 				sequence_dict = {'header': line[1:], 'sequence': ''}
# 			else:
# 				sequence_dict['sequence'] += line
#
# 	return sequence_dict
#
#
# # NOT_USED
# def compute_consensus_sequence(reference_file, sequence_to_analyse, compressed_vcf_file, outdir, sufix):
# 	sequence_dict = None
#
# 	gene_fasta = os.path.join(outdir, str(sequence_to_analyse + '.fasta'))
#
# 	run_successfully, stdout = index_fasta_samtools(reference_file, sequence_to_analyse, gene_fasta, False)
# 	if run_successfully:
# 		command = ['bcftools', 'consensus', '-f', gene_fasta, compressed_vcf_file]
# 		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
# 		if run_successfully:
# 			sequence_dict = parse_fasta_inMemory(stdout)
#
# 	return run_successfully, sequence_dict


def create_sample_consensus_sequence(outdir, sequence_to_analyse, reference_file, variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence, length_extra_seq):
	variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found = get_true_variants(variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence)

	variants_correct, variants_noMatter, variants_alignment, number_multi_alleles = cleanning_variants_extra_seq(variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found, length_extra_seq, len(sequence))

	# run_successfully = False
	# consensus = {'correct': {}, 'noMatter': {}, 'alignment': {}}
	# for variant_type in ['variants_correct', 'variants_noMatter', 'variants_alignment']:
	# 	run_successfully, compressed_vcf_file = write_variants_vcf(eval(variant_type), outdir, sequence_to_analyse, variant_type.split('_', 1)[1])
	# 	if run_successfully:
	# 		run_successfully, sequence_dict = compute_consensus_sequence(reference_file, sequence_to_analyse, compressed_vcf_file, outdir, variant_type.split('_', 1)[1])
	# 		if run_successfully:
	# 			consensus[variant_type.split('_', 1)[1]] = {'header': sequence_dict['header'], 'sequence': sequence_dict['sequence'][length_extra_seq:len(sequence_dict['sequence']) - length_extra_seq]}
	#
	# return run_successfully, number_multi_alleles, consensus
	return True, number_multi_alleles, None


@utils.trace_unhandled_exceptions
def analyse_sequence_data(bam_file, sequence_information, outdir, counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	percentage_absent = None
	percentage_lowCoverage = None
	meanCoverage = None

	# Create vcf file (for multiple alleles check)
	run_successfully, gene_vcf = create_vcf(bam_file, sequence_information['header'], outdir, counter, reference_file)

	if run_successfully:
		# Create coverage tab file
		run_successfully, gene_coverage = compute_genome_coverage_data(bam_file, sequence_information['header'], outdir, counter)

		if run_successfully:
			variants = get_variants(gene_vcf)

			coverage = get_coverage(gene_coverage)

			run_successfully, number_multi_alleles, consensus_sequence = create_sample_consensus_sequence(outdir, sequence_information['header'], reference_file, variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence_information['sequence'], length_extra_seq)

			percentage_absent, percentage_lowCoverage, meanCoverage = get_coverage_report(coverage, sequence_information['length'], minimum_depth_presence, minimum_depth_call, length_extra_seq)

	utils.saveVariableToPickle([run_successfully, counter, number_multi_alleles, percentage_absent, percentage_lowCoverage, meanCoverage, consensus_sequence], outdir, str('coverage_info.' + str(counter)))


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
						sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': '', 'length': 0}
					else:
						sequence_dict[sequence_counter]['sequence'] += line
						sequence_dict[sequence_counter]['length'] += len(line)
				else:
					sequence_dict = None
			else:
				blank_line_found = True

	return sequence_dict


def sequence_data(sample, reference_file, bam_file, outdir, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	sequence_data_outdir = os.path.join(outdir, 'sequence_data', '')
	utils.removeDirectory(sequence_data_outdir)
	os.mkdir(sequence_data_outdir)

	sequences = get_sequence_information(reference_file)

	pool = multiprocessing.Pool(processes=threads)
	for sequence_counter in sequences:
		sequence_dir = os.path.join(sequence_data_outdir, str(sequence_counter), '')
		utils.removeDirectory(sequence_dir)
		os.makedirs(sequence_dir)
		pool.apply_async(analyse_sequence_data, args=(bam_file, sequences[sequence_counter], sequence_dir, sequence_counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele,))
	pool.close()
	pool.join()

	run_successfully, sample_data = gather_data_together(sample, sequence_data_outdir, sequences, outdir.rsplit('/', 2)[0])

	return run_successfully, sample_data


def chunkstring(string, length):
	return (string[0 + i:length + i] for i in range(0, len(string), length))


# # NOT_USED
# def write_consensus(outdir, sample, consensus_sequence):
# 	for consensus_type in ['correct', 'noMatter', 'alignment']:
# 		with open(os.path.join(outdir, str(sample + '.' + consensus_type + '.fasta')), 'at') as writer:
# 			writer.write('>' + consensus_sequence[consensus_type]['header'] + '\n')
# 			fasta_sequence_lines = chunkstring(consensus_sequence[consensus_type]['sequence'], 80)
# 			for line in fasta_sequence_lines:
# 				writer.write(line + '\n')


def gather_data_together(sample, data_directory, sequences_information, outdir):
	run_successfully = True
	counter = 0
	sample_data = {}

	# write_consensus_first_time = True

	genes_directories = [d for d in os.listdir(data_directory) if not d.startswith('.') and os.path.isdir(os.path.join(data_directory, d, ''))]
	for gene_dir in genes_directories:
		gene_dir_path = os.path.join(data_directory, gene_dir, '')

		files = [f for f in os.listdir(gene_dir_path) if not f.startswith('.') and os.path.isfile(os.path.join(gene_dir_path, f))]
		for file_found in files:
			if file_found.startswith('coverage_info.') and file_found.endswith('.pkl'):
				file_path = os.path.join(gene_dir_path, file_found)

				if run_successfully:
					run_successfully, sequence_counter, multiple_alleles_found, percentage_absent, percentage_lowCoverage, meanCoverage, consensus_sequence = utils.extractVariableFromPickle(file_path)

					# if write_consensus_first_time:
					# 	for consensus_type in ['correct', 'noMatter', 'alignment']:
					# 		file_to_remove = os.path.join(outdir, str(sample + '.' + consensus_type + '.fasta'))
					# 		if os.path.isfile(file_to_remove):
					# 			os.remove(file_to_remove)
					# 	write_consensus_first_time = False
					# write_consensus(outdir, sample, consensus_sequence)

					sample_data[sequence_counter] = {'header': sequences_information[sequence_counter]['header'], 'gene_coverage': 100 - percentage_absent, 'gene_low_coverage': percentage_lowCoverage, 'gene_number_positions_multiple_alleles': multiple_alleles_found, 'gene_mean_read_coverage': meanCoverage}
					counter += 1

		utils.removeDirectory(gene_dir_path)

	if counter != len(sequences_information):
		run_successfully = False

	return run_successfully, sample_data


# rematch_timer = functools.partial(utils.timer, name='ReMatCh module')
trueCoverage_timer = functools.partial(utils.timer, name='True coverage check')


# @rematch_timer
# def runTrueCoverage(sample, fastq_files, reference_file, threads, outdir, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, minimum_gene_coverage, conserved_True):
@trueCoverage_timer
def runTrueCoverage(sample, fastq_files, reference_file, threads, outdir, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, minimum_gene_coverage, maximum_number_absent_genes, maximum_number_genes_multiple_alleles, minimum_read_coverage, conserved_True):
	pass_qc = False
	failing = {}

	trueCoverage_folder = os.path.join(outdir, 'trueCoverage', '')
	utils.removeDirectory(trueCoverage_folder)
	os.mkdir(trueCoverage_folder)

	# Map reads
	run_successfully, bam_file, reference_file = mapping_reads(fastq_files, reference_file, threads, trueCoverage_folder, conserved_True)

	if run_successfully:
		# Index reference file
		run_successfully, stdout = index_fasta_samtools(reference_file, None, None, True)
		if run_successfully:
			print 'Analysing alignment data'
			run_successfully, sample_data = sequence_data(sample, reference_file, bam_file, trueCoverage_folder, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

			if run_successfully:
				print 'Writing report file'
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

					if len(sample_data) - number_absent_genes > 0:
						mean_sample_coverage = float(mean_sample_coverage) / float(len(sample_data) - number_absent_genes)
					else:
						mean_sample_coverage = 0

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

	# return run_successfully, sample_data if 'sample_data' in locals() else None, {'number_absent_genes': number_absent_genes, 'number_genes_multiple_alleles': number_genes_multiple_alleles, 'mean_sample_coverage': round(mean_sample_coverage, 2)} if 'number_absent_genes' in locals() else None

	return run_successfully, pass_qc, failing
