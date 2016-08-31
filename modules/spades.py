import utils
import os
import shutil
from functools import partial


# Run Spades
def spades(spades_folder, threads, fastq_files, notUseCareful, maxMemory, minCoverageAssembly, kmers):
	contigs = os.path.join(spades_folder, 'contigs.fasta')

	command = ['spades.py', '', '--only-assembler', '--threads', str(threads), '--memory', str(maxMemory), '--cov-cutoff', str(minCoverageAssembly), '', '-1', fastq_files[0], '-2', fastq_files[1], '-o', spades_folder]

	if not notUseCareful:
		command[1] = '--careful'

	if len(kmers) > 0:
		kmers = ','.join(map(str, kmers))
		command[9] = str('-k ' + kmers)

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	return run_successfully, contigs


# Rename contigs and contigs.fasta file while filtering for contigs length
def renameFilterContigs(sampleName, outdir, spadesContigs, minContigsLength, minCoverageContigs):
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
				items = contigHeader.split('_')
				if len(contigSequence) >= minContigsLength and float(items[5]) >= minCoverageContigs:
					writer.write(">" + sampleName + "_" + contigHeader + "\n")
					writer.write(contigSequence + "\n")
					number_bases = number_bases + len(contigSequence)
					number_contigs = number_contigs + 1
			contigHeader = ""
			contigSequence = ""
			contigHeader = line[1:].splitlines()[0]
		else:
			contigSequence = contigSequence + line.splitlines()[0]
	items = contigHeader.split('_')
	if len(contigSequence) >= minContigsLength and float(items[5]) >= minCoverageContigs:
		writer.write(">" + sampleName + "_" + contigHeader + "\n")
		writer.write(contigSequence + "\n")
		number_bases = number_bases + len(contigSequence)
		number_contigs = number_contigs + 1
	writer.close()

	return newContigsFile, number_contigs, number_bases


def define_kmers(kmers, maximumReadsLength):
	kmers_use = []
	if maximumReadsLength is not None:
		for kmer in kmers:
			if kmer <= maximumReadsLength:
				kmers_use.append(kmer)
	return kmers_use


spades_timer = partial(utils.timer, name='SPAdes')


# Run SPAdes procedure
@spades_timer
def runSpades(sampleName, outdir, threads, fastq_files, notUseCareful, maxMemory, minCoverageAssembly, minContigsLength, estimatedGenomeSizeMb, kmers, maximumReadsLength, saveReport, defaultKmers, minCoverageContigs):
	pass_qc = False
	failing = {}
	failing['sample'] = False

	contigs = None

	# Create SPAdes output directory
	spades_folder = os.path.join(outdir, 'spades', '')
	utils.removeDirectory(spades_folder)
	os.mkdir(spades_folder)

	if defaultKmers:
		kmers = []
	else:
		kmers = define_kmers(kmers, maximumReadsLength)
		if len(kmers) == 0:
			print 'SPAdes will use its default k-mers'
		else:
			print 'SPAdes will use the following k-mers: ' + str(kmers)

	run_successfully, contigs = spades(spades_folder, threads, fastq_files, notUseCareful, maxMemory, minCoverageAssembly, kmers)

	if run_successfully:
		shutil.copyfile(contigs, os.path.join(outdir, 'SPAdes_original_assembly.contigs.fasta'))
		print 'Filtering for contigs with at least ' + str(minContigsLength) + ' nucleotides and a coverage of ' + str(minCoverageContigs)
		contigsFiltered, number_contigs, number_bases = renameFilterContigs(sampleName, outdir, contigs, minContigsLength, minCoverageContigs)
		print str(number_bases) + ' assembled nucleotides in ' + str(number_contigs) + ' contigs'

		if saveReport:
			report_file = os.path.join(outdir, 'spades_report.txt')
			with open(report_file, 'wt') as writer:
				writer.write('#contigs' + '\n' + str(number_contigs) + '\n' + '#bp' + '\n' + str(number_bases) + '\n')
				writer.flush()

		if number_bases >= estimatedGenomeSizeMb * 1000000 * 0.8 and number_bases <= estimatedGenomeSizeMb * 1000000 * 1.5:
			if number_contigs <= 100 * number_bases / 1500000:
				pass_qc = True
			else:
				failing['sample'] = 'The number of assembled contigs (' + str(number_contigs) + ') exceeds ' + str(100 * number_bases / 1500000)
				print failing['sample']
		else:
			failing['sample'] = 'The number of assembled nucleotides (' + str(number_bases) + ') are lower than 80% or higher than 150% of the provided estimated genome size'
			print failing['sample']
	else:
		failing['sample'] = 'Did not run'
		print failing['sample']
		contigsFiltered = None

	utils.removeDirectory(spades_folder)

	return run_successfully, pass_qc, failing, contigsFiltered
