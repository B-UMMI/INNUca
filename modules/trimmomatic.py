import os
import utils
from functools import partial


# Run Trimmomatic
def trimmomatic(sampleName, trimmomatic_folder, threads, adaptersFasta, script_path, doNotSearchAdapters, fastq_files, maxReadsLength, doNotTrimCrops, crop, headCrop, leading, trailing, slidingWindow, minLength, nts2clip_based_ntsContent):
	fastq = sorted(fastq_files)[0]

	# Run Trimmomatic
	command = ['trimmomatic-0.36.jar', 'PE', '-threads', str(threads), '-basein', fastq, '-baseout', os.path.join(trimmomatic_folder, str(sampleName + '.fastq.gz')), '', '', '', str('SLIDINGWINDOW:' + slidingWindow), str('LEADING:' + str(leading)), str('TRAILING:' + str(trailing)), str('MINLEN:' + str(minLength)), 'TOPHRED33', '']

	if not doNotTrimCrops:
		if maxReadsLength is not None:
			if crop is not None:
				crop = maxReadsLength - crop[0]
				command[8] = str('CROP:' + str(crop))
			else:
				if nts2clip_based_ntsContent is not None:
					crop = nts2clip_based_ntsContent[1]
					print str(crop) + ' nucleotides will be clipped at the end of reads'
					crop = maxReadsLength - crop
					command[8] = str('CROP:' + str(crop))
		else:
			print 'Because FastQC did not run successfully, --trimCrop option will not be considered'

		if headCrop is not None:
			command[9] = str('HEADCROP:' + str(headCrop[0]))
		else:
			if nts2clip_based_ntsContent is not None:
				headCrop = nts2clip_based_ntsContent[0]
				print str(headCrop) + ' nucleotides will be clipped at the beginning of reads'
				command[9] = str('HEADCROP:' + str(headCrop))

	if not doNotSearchAdapters:
		if adaptersFasta is not None:
			print 'Removing adapters contamination using ' + adaptersFasta
			command[10] = 'ILLUMINACLIP:' + adaptersFasta + ':3:30:10:6:true'
		else:
			trimmomatic_adapters_folder = os.path.join(os.path.dirname(script_path), 'src', 'Trimmomatic-0.36', 'adapters')
			adapters_files = [os.path.join(trimmomatic_adapters_folder, 'NexteraPE-PE.fa'), os.path.join(trimmomatic_adapters_folder, 'TruSeq2-PE.fa'), os.path.join(trimmomatic_adapters_folder, 'TruSeq3-PE-2.fa')]
			print 'Removing adapters contamination using ' + str(adapters_files)
			adaptersFasta = concatenateFastaFiles(adapters_files, trimmomatic_folder, 'concatenated_adaptersFile.fasta')
			command[10] = 'ILLUMINACLIP:' + adaptersFasta + ':3:30:10:6:true'

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	if not run_successfully:
		print 'Trimmomatic fail! Trying run with Phred+33 enconding defined...'
		command[16] = '-phred33'
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)
		if not run_successfully:
			print 'Trimmomatic fail again! Trying run with Phred+64 enconding defined...'
			command[16] = '-phred64'
			run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	return run_successfully


# Concatenate fasta files
def concatenateFastaFiles(list_fasta_files, outdir, outFileName):
	concatenated_file = os.path.join(outdir, outFileName)
	with open(concatenated_file, 'wt') as outfile:
		for fasta in list_fasta_files:
			with open(fasta, 'rtU') as infile:
				for line in infile:
					if len(line) > 0:
						line = line.splitlines()[0]
						outfile.write(line + '\n')
						outfile.flush()

	return concatenated_file


# Scans the trimmed files for paired and unpaired reads
def getTrimmomaticPairedReads(trimmomatic_folder):
	paired_reads = []

	files = [f for f in os.listdir(trimmomatic_folder) if not f.startswith('.') and os.path.isfile(os.path.join(trimmomatic_folder, f))]

	for fastq in files:
		if fastq.endswith('P.fastq.gz'):
			paired_reads.append(os.path.join(trimmomatic_folder, fastq))

	return paired_reads


# Chech whether Trimmomatic returns 0 reads
def controlForZeroReads(fastq_files):
	not_empty_fastq = False

	fastq = fastq_files[0]

	command = ['', '--stdout', '--keep', fastq, '|', 'head', '-n', '4']

	filetype = utils.compressionType(fastq)
	if filetype == 'gz':
		command[0] = 'gunzip'
	elif filetype == 'bz2':
		command[0] = 'bunzip2'

	if command[0] != '':
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None)

		if run_successfully:
			stdout = stdout.splitlines()
			if len(stdout) == 4:
				not_empty_fastq = True

	return not_empty_fastq


trim_timer = partial(utils.timer, name='Trimmomatic')


# Run Trimmomatic procedure
@trim_timer
def runTrimmomatic(sampleName, outdir, threads, adaptersFasta, script_path, doNotSearchAdapters, fastq_files, maxReadsLength, doNotTrimCrops, crop, headCrop, leading, trailing, slidingWindow, minLength, nts2clip_based_ntsContent):
	failing = {}
	failing['sample'] = False
	not_empty_fastq = False

	paired_reads = None

	# Create Trimmomatic output directory
	trimmomatic_folder = os.path.join(outdir, 'trimmomatic', '')
	utils.removeDirectory(trimmomatic_folder)
	os.mkdir(trimmomatic_folder)

	run_successfully = trimmomatic(sampleName, trimmomatic_folder, threads, adaptersFasta, script_path, doNotSearchAdapters, fastq_files, maxReadsLength, doNotTrimCrops, crop, headCrop, leading, trailing, slidingWindow, minLength, nts2clip_based_ntsContent)

	if run_successfully:
		paired_reads = getTrimmomaticPairedReads(trimmomatic_folder)
		not_empty_fastq = controlForZeroReads(paired_reads)

		if not not_empty_fastq:
			failing['sample'] = 'Zero reads after Trimmomatic'
			print failing['sample']
	else:
		failing['sample'] = 'Did not run'
		print failing['sample']

	return run_successfully, not_empty_fastq, failing, paired_reads, trimmomatic_folder
