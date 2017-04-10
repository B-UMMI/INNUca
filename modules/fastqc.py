import sys
import os
import utils
from functools import partial


# Prepare adapters fasta file to FastQC (tabular file)
def adapters2fastQC(outdir, adaptersFastaFile):
	adaptersFile = os.path.join(outdir, 'temp.adapters2fastQC.tab')
	writer = open(adaptersFile, 'wt')
	adapterHeader = ''
	adapterSequence = ''
	with open(adaptersFastaFile, 'rtU') as adapters:
		for line in adapters:
			if line.startswith('>'):
				if adapterHeader != '':
					writer.write(adapterHeader + '\t' + adapterSequence + '\n')
				adapterHeader = ''
				adapterSequence = ''
				adapterHeader = line[1:].splitlines()[0]
			else:
				adapterSequence = adapterSequence + line.splitlines()[0]
		writer.write(adapterHeader + '\t' + adapterSequence + '\n')
	writer.close()
	return adaptersFile


# Run FastQC
def fastQC(fastqc_folder, threads, adaptersFasta, fastq_files):
	# Create temporary FastQC foldes
	os.mkdir(os.path.join(fastqc_folder, 'temp.fastqc_temporary_dir', ''))

	# Run FastQC
	command = ['fastqc', '-o', fastqc_folder, '--extract', '--nogroup', '--format', 'fastq', '--threads', str(threads), '', '--dir', os.path.join(fastqc_folder, 'temp.fastqc_temporary_dir', '')]
	command = command + fastq_files
	if adaptersFasta is not None:
		adaptersTEMP = adapters2fastQC(fastqc_folder, adaptersFasta)
		print 'Scanning for adapters contamination using ' + adaptersFasta
		command[9] = '--adapters ' + adaptersTEMP
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

	# Remove temporary files
	os.rmdir(os.path.join(fastqc_folder, 'temp.fastqc_temporary_dir', ''))
	if adaptersFasta is not None:
		os.remove(adaptersTEMP)

	return run_successfully


# Parse FastQC run
def parseFastQC(fastqc_folder, fastq_files):
	fastq_files = map(os.path.basename, fastq_files)
	goodReads = []
	badReads = []
	failing = {}
	warning = {}
	for reads in fastq_files:
		reads_file = reads.rsplit('.', 2)[0]
		with open(os.path.join(fastqc_folder, str(reads_file + '_fastqc'), 'summary.txt'), 'rtU') as fastqc_summary:
			fastqc = {}
			for line in fastqc_summary:
				fastqc[line.split('\t')[1]] = line.split()[0]

			bad_fastq = False
			failing[reads] = []

			if fastqc.get('Per base sequence quality') != 'PASS':
				bad_fastq = True
				failing[reads].append('Bad per base sequence quality')
			if fastqc.get('Per sequence GC content') == 'FAIL':
				bad_fastq = True
				failing[reads].append('Bad per sequence GC content')
			if fastqc.get('Per base N content') != 'PASS':
				bad_fastq = True
				failing[reads].append('Bad per base N content')
			if fastqc.get('Sequence Length Distribution') == 'FAIL':
				bad_fastq = True
				failing[reads].append('Bad sequence length distribution')
			if fastqc.get('Overrepresented sequences') != 'PASS':
				bad_fastq = True
				failing[reads].append('Overrepresented sequences')
			if fastqc.get('Adapter Content') != 'PASS':
				bad_fastq = True
				failing[reads].append('Found adapters sequences')

			if fastqc.get('Per base sequence content') == 'FAIL':
				if reads not in warning:
					warning[reads] = []
				warning[reads].append('Bad per base sequence content')

			if not bad_fastq:
				goodReads.append(reads)
			else:
				badReads.append(reads)

	for fastq in failing:
		if len(failing[fastq]) == 0:
			failing[fastq] = False

	return goodReads, badReads, failing, warning


# Get reads length data, nucleotide bias status & number of reads
def getReadsInformation(fastqc_folder, fastq_files):
	fastq_files = map(os.path.basename, fastq_files)
	maxReadsLengths = []
	moreFrequentReadsLength = []
	numberReads = []
	ntsContent_biasStatus = {}
	for reads in fastq_files:
		reads_file = reads.rsplit('.', 2)[0]
		with open(os.path.join(fastqc_folder, str(reads_file + '_fastqc'), 'fastqc_data.txt'), 'rtU') as fastqc_data:
			modulesAssessed = 0

			basicModule = False
			maxLength = None

			nt_content_Module = False
			ntsContent_biasStatus[reads] = []

			lengthModule = False
			numberReadsPerLength = None
			readsLength = None
			for line in fastqc_data:
				if not line.startswith('#'):
					if line.startswith('>>'):
						if line.startswith('Basic Statistics', 2):
							basicModule = True
						elif line.startswith('Per base sequence content', 2):
							nt_content_Module = True
						elif line.startswith('Sequence Length Distribution', 2):
							lengthModule = True
						else:
							if basicModule:
								basicModule = False
								modulesAssessed += 1
							elif nt_content_Module:
								nt_content_Module = False
								modulesAssessed += 1
							elif lengthModule:
								lengthModule = False
								modulesAssessed += 1
							else:
								if modulesAssessed == 3:
									break
					else:
						if basicModule:
							line_splited = line.splitlines()[0].split('\t')
							if line_splited[0] == 'Sequence length':
								lengths = line_splited[1].split('-')
								if len(lengths) == 1:
									maxLength = int(lengths[0])
								elif len(lengths) == 2:
									maxLength = int(lengths[1])
								else:
									sys.exit("ERROR: strange FastQC 'Sequence length' information")
							elif line_splited[0] == 'Total Sequences':
								numberReads.append(int(line_splited[1]))
						# Get nucleotide content bias status for read position
						elif nt_content_Module:
							line_splited = line.splitlines()[0].split('\t')
							gc = (float(line_splited[1]) + 0.1) / (float(line_splited[4]) + 0.1)
							at = (float(line_splited[2]) + 0.1) / (float(line_splited[3]) + 0.1)
							if 0.8 <= gc <= 1.2 and 0.8 <= at <= 1.2:
								ntsContent_biasStatus[reads].append('unbiased')
							else:
								ntsContent_biasStatus[reads].append('biased')
						elif lengthModule:
							line_splited = line.splitlines()[0].split('\t')
							number_reads = float(line_splited[1])
							if numberReadsPerLength < number_reads:
								numberReadsPerLength = number_reads
								readsLength = line_splited[0]

			maxReadsLengths.append(maxLength)

			moreFrequentReadsLength.append(readsLength)

	maximumReadsLength = max(maxReadsLengths)

	# maximumReadsLength for both fastq files
	# moreFrequentReadsLength for each fastq file in the format of 50-54
	# numberReads for each fastq file
	return maximumReadsLength, moreFrequentReadsLength, numberReads, ntsContent_biasStatus


# Get the number of nucleotides for 5' and 3' that we want trimmomatic to clip based on nucleotide content bias for each fastq file
def nts2clip(dict_fastqs_ntsBiased_status):
	nts2clip_based_ntsContent = {}

	for fastq in dict_fastqs_ntsBiased_status:
		nts2clip_based_ntsContent[fastq] = [0, 0]
		nt_content = dict_fastqs_ntsBiased_status[fastq]
		five_end = False
		reads_range = list(range(0, len(nt_content) - 1))
		for i in reads_range:
			if nt_content[i] == 'biased':
				if i <= len(nt_content) / 2:
					if not five_end:
						if nt_content[i + 1] == 'unbiased' and nt_content[i + 2] == 'unbiased':
							nts2clip_based_ntsContent[fastq][0] = i + 1
							five_end = True
				else:
					if i + 1 not in reads_range:
							nts2clip_based_ntsContent[fastq][1] = len(nt_content) - i
							break
					elif i + 2 not in reads_range:
							nts2clip_based_ntsContent[fastq][1] = len(nt_content) - i
							break
					else:
						if nt_content[i + 1] == 'unbiased' and nt_content[i + 2] == 'unbiased':
							continue
						else:
							nts2clip_based_ntsContent[fastq][1] = len(nt_content) - i
							break

	nts2clip_based_ntsContent = [max(nts2clip_based_ntsContent[nts2clip_based_ntsContent.keys()[0]][0], nts2clip_based_ntsContent[nts2clip_based_ntsContent.keys()[1]][0]), min(nts2clip_based_ntsContent[nts2clip_based_ntsContent.keys()[0]][1], nts2clip_based_ntsContent[nts2clip_based_ntsContent.keys()[1]][1])]

	return nts2clip_based_ntsContent


def check_FastQC_runSuccessfully(fastqc_folder, fastq_files):
	run_successfully = True
	for reads in fastq_files:
		reads_file = os.path.basename(reads).rsplit('.', 2)[0]
		try:
			open(os.path.join(fastqc_folder, str(reads_file + '_fastqc'), 'fastqc_data.txt'), 'rtU')
		except Exception as e:
			print e
			run_successfully = False
	return run_successfully


fastqc_timer = partial(utils.timer, name='FastQC analysis')


# Run FastQC analysis
@fastqc_timer
def runFastQCanalysis(outdir, threads, adaptersFasta, fastq_files, keepFiles, fastQC_run_name):
	pass_qc = False
	failing = {}
	failing['sample'] = False

	warnings = {}

	maximumReadsLength = None
	nts2clip_based_ntsContent = None

	# Create FastQC output directory
	fastqc_folder = os.path.join(outdir, str('fastqc_' + fastQC_run_name), '')
	utils.removeDirectory(fastqc_folder)
	os.mkdir(fastqc_folder)

	# Run FastQC
	run_successfully = fastQC(fastqc_folder, threads, adaptersFasta, fastq_files)
	if run_successfully:
		# Check whether FastQC really run_successfully
		run_successfully = check_FastQC_runSuccessfully(fastqc_folder, fastq_files)
		if not run_successfully:
			failing['sample'] = 'Did not run'
			return run_successfully, pass_qc, failing, warnings, maximumReadsLength, nts2clip_based_ntsContent

		# Check which reads pass FastQC
		goodReads, badReads, failing, warnings = parseFastQC(fastqc_folder, fastq_files)
		# Get reads information
		maximumReadsLength, moreFrequentReadsLength, numberReads, ntsContent_biasStatus = getReadsInformation(fastqc_folder, fastq_files)
		# Get number nucleotides to clip based on nucleotide content bias
		nts2clip_based_ntsContent = nts2clip(ntsContent_biasStatus)

		print "Number of reads found: " + str(numberReads)
		print "Maximum reads length found for both fastq files: " + str(maximumReadsLength) + " nts"
		print "Reads length class more frequently found in fastq files: " + str(moreFrequentReadsLength)
		if len(badReads) == 0:
			pass_qc = True
		elif len(badReads) > 0:
			print "Reads files FAILING FastQC control: " + str(badReads)
		if len(goodReads) > 0:
			print "Reads files passing FastQC control: " + str(goodReads)
		print 'To improve reads quality, consider clipping the next number of nucleotides in the fastq files at 5 end and 3 end, respectively: ' + str(nts2clip_based_ntsContent)
	else:
		failing['sample'] = 'Did not run'
		print failing['sample']

	if not keepFiles:
		utils.removeDirectory(fastqc_folder)

	return run_successfully, pass_qc, failing, warnings, maximumReadsLength, nts2clip_based_ntsContent
