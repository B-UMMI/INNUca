import utils
import os
import shutil


# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)
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
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	if not run_successfully:
		sam_file = None

	return run_successfully, sam_file


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
	outFormat_string = os.path.splitext(output_file)[1][1:].lower()
	command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
	if sortByName_True:
		command[6] = '-n'
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)
	if not run_successfully:
		output_file = None
	return run_successfully, output_file


# Index alignment file
def indexAlignment(alignment_file):
	command = ['samtools', 'index', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)
	return run_successfully


def pilon(assembly, bam_file, outdir):
	assembly_polished = os.path.splitext(assembly)[0] + '.polished.fasta'
	command = ['pilon-1.18.jar', '--genome', assembly, '--frags', bam_file, '--outdir', outdir, '--output', os.path.basename(os.path.splitext(assembly_polished)[0]), '--changes', '--vcf']
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)
	if not run_successfully:
		assembly_polished = None
	return run_successfully, assembly_polished


def parsePilonResult(assembly_polished, outdir):
	corrections = {}
	with open(str(os.path.splitext(assembly_polished)[0] + '.changes'), 'rtU') as reader:
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				contig = line.split(' ', 1)[0].rsplit(':', 1)[0]
				if contig in corrections:
					corrections[contig] += 1
				else:
					corrections[contig] = 1
	report_file = os.path.join(outdir, 'pilon_report.txt')
	with open(report_file, 'wt') as writer:
		writer.write('#general' + '\n')
		writer.write('>changes' + '\n' + str(sum(corrections.values())) + '\n' + '>contigs' + '\n' + str(len(corrections)) + '\n')
		writer.flush()
		writer.write('#by_contigs' + '\n')
		for contig in corrections:
			writer.write(str('>' + contig) + '\n' + str(corrections[contig]) + '\n')
	print str(sum(corrections.values())) + ' changes made by Pilon in ' + str(len(corrections)) + ' contigs'


# Count sequenced bases
def runPilon(assembly, fastq_files, threads, outdir, keepFiles, keepSPAdesAssembly):
	failing = {}
	failing['sample'] = False

	pilon_folder = os.path.join(outdir, 'pilon', '')
	utils.removeDirectory(pilon_folder)
	os.mkdir(pilon_folder)

	# Create a symbolic link to the assembly
	assembly_link = os.path.join(pilon_folder, os.path.basename(assembly))
	os.symlink(assembly, assembly_link)

	assembly_polished = None

	# Index assembly using Bowtie2
	run_successfully = indexSequenceBowtie2(assembly_link, threads)

	if run_successfully:
		run_successfully, sam_file = mappingBowtie2(fastq_files, assembly_link, threads, pilon_folder)

		if run_successfully:
			bam_file = os.path.splitext(sam_file)[0] + '.bam'
			run_successfully, bam_file = sortAlignment(sam_file, bam_file, False, threads)

			if run_successfully:
				os.remove(sam_file)
				run_successfully = indexAlignment(bam_file)

				if run_successfully:
					run_successfully, assembly_polished = pilon(assembly_link, bam_file, pilon_folder)

					if run_successfully:
						parsePilonResult(assembly_polished, outdir)
						shutil.copyfile(assembly_polished, os.path.join(outdir, os.path.basename(assembly_polished)))
						assembly_polished = os.path.join(outdir, os.path.basename(assembly_polished))
						if not keepSPAdesAssembly:
							os.remove(assembly)

	if not keepFiles:
		utils.removeDirectory(pilon_folder)

	return run_successfully, failing, assembly_polished
