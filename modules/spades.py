import utils
import os
import shutil

# Run Spades
def spades(spades_folder, threads, fastq_files, notUseCareful, maxMemory, minCoverage, kmers):
	contigs = os.path.join(spades_folder, 'contigs.fasta')

	command = ['spades.py', '', '--only-assembler', '--threads', str(threads), '--memory', str(maxMemory), '--cov-cutoff', str(minCoverage), '', '-1', fastq_files[0], '-2', fastq_files[1], '-o', spades_folder]

	if not notUseCareful:
		command[1] = '--careful'

	if len(kmers) > 0:
		kmers = ','.join(map(str, kmers))
		command[9] = str('-k ' + kmers)

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command)

	return run_successfully, contigs

# Rename contigs and contigs.fasta file while filtering for contigs length
def renameFilterContigs(sampleName, outdir, spadesContigs, minContigsLength):
	newContigsFile = os.path.join(outdir, str(sampleName + '.contigs.fasta'))
	number_contigs = 0
	number_bases = 0

	writer = open(newContigsFile, 'wt')
	contigHeader = ""
	contigSequence = ""
	contigs = open(spadesContigs)
	for line in contigs:
		if line[0] == '>':
			if contigHeader != "":
				if len(contigSequence) >= minContigsLength:
					writer.write(contigHeader + "\n")
					writer.write(contigSequence + "\n")
					number_bases = number_bases + len(contigSequence)
					number_contigs = number_contigs + 1
			contigHeader = ""
			contigSequence = ""
			contigHeader = ">" + sampleName + "_" + line[1:].splitlines()[0]
		else:
			contigSequence = contigSequence + line.splitlines()[0]
	if len(contigSequence) >= minContigsLength:
		writer.write(contigHeader + "\n")
		writer.write(contigSequence + "\n")
		number_bases = number_bases + len(contigSequence)
		number_contigs = number_contigs + 1
	writer.close()

	return newContigsFile, number_contigs, number_bases

def define_kmers(kmers, maximumReadsLength):
	if isinstance(kmers[0], (list, tuple)):
		kmers = kmers[0]
	kmers_use = []
	if maximumReadsLength is not None:
		for kmer in kmers:
			if kmer <= maximumReadsLength:
				kmers_use.append(kmer)
	return kmers_use

# Run SPAdes procedure
def runSpades(sampleName, outdir, threads, fastq_files, notUseCareful, maxMemory, minCoverage, minContigsLength, estimatedGenomeSizeMb, kmers, maximumReadsLength, saveReport):
	pass_qc = False
	failing = {}
	failing['sample'] = False

	contigs = None

	# Create SPAdes output directory
	spades_folder = os.path.join(outdir, 'spades', '')
	utils.removeDirectory(spades_folder)
	os.mkdir(spades_folder)

	kmers = define_kmers(kmers, maximumReadsLength)
	if len(kmers) == 0:
		print 'SPAdes will use its default k-mers'
	else:
		print 'SPAdes will use the following k-mers: ' + str(kmers)

	run_successfully, contigs = spades(spades_folder, threads, fastq_files, notUseCareful, maxMemory, minCoverage, kmers)

	if run_successfully:
		print 'Filtering for contigs with at least ' + str(minContigsLength) + ' nucleotides'
		contigsFiltered, number_contigs, number_bases = renameFilterContigs(sampleName, outdir, contigs, minContigsLength)
		print str(number_bases) + ' assembled nucleotides in ' + str(number_contigs) + ' contigs'

		if saveReport:
			report_file = os.path.join(outdir, 'spades_report.txt')
			with open(report_file, 'wt') as writer:
				writer.write('#contigs' + '\n' + str(number_contigs) + '\n' + '#bp' + '\n' + str(number_bases) + '\n')
				writer.flush()

		if number_bases >= estimatedGenomeSizeMb*1000000*0.8 and number_bases <= estimatedGenomeSizeMb*1000000*1.5:
			if number_contigs <= 100*number_bases/1500000:
				pass_qc = True
			else:
				failing['sample'] = 'The number of assembled contigs (' + str(number_contigs) + ') exceeds ' + str(100*number_bases/1500000)
				print failing['sample']
		else:
			failing['sample'] = 'The number of assembled nucleotides (' + str(number_bases) + ') are lower than 80% or higher than 150% of the provided estimated genome size'
			print failing['sample']
	else:
		failing['sample'] = 'Did not run'
		print failing['sample']

	utils.removeDirectory(spades_folder)

	return run_successfully, pass_qc, failing, contigsFiltered
