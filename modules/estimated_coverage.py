import utils
import os
from functools import partial

coverage_timer = partial(utils.timer, name='First Estimated Coverage analysis')

# Count sequenced bases
def countSequencedBases(fastq_file):
	run_successfully = False
	bases = None

	command = ['', '--keep', '--stdout', fastq_file, '|', 'grep', '--after-context=1', '"@"', '|', 'grep', '--invert-match', '"^--$"', '|', 'grep', '--invert-match', '"@"', '|', 'wc', '']

	# Determine compression type
	filetype = utils.compressionType(fastq_file)
	if filetype == 'gz':
		command[0] = 'gunzip'
	elif filetype == 'bz2':
		command[0] = 'bunzip2'
	else:
		return run_successfully, bases

	# Number of characters
	command[18] = '--chars'
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None)
	if run_successfully:
		bases = int(stdout.splitlines()[0])

		# Number of lines
		command[18] = '--lines'
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None)
		if run_successfully:
			lines = int(stdout.splitlines()[0])
			bases = bases - lines

	return run_successfully, bases


# Get extimated coverage

@coverage_timer
def getEstimatedCoverage(fastq_files, estimatedGenomeSizeMb, outdir):
	run_successfully = False
	pass_qc = False
	failing = {}
	failing['sample'] = False

	# Run Estimated Coverage
	numberBases = 0
	for fastq in fastq_files:
		# Get number bases for each fastq file
		run_successfully, bases = countSequencedBases(fastq)
		if not run_successfully:
			break
		else:
			numberBases = numberBases + bases

	if run_successfully:
		estimatedCoverage = numberBases / float(estimatedGenomeSizeMb * 1000000)
		estimatedCoverage = round(estimatedCoverage, 1)

		report_file = os.path.join(outdir, 'coverage_report.txt')
		report = str(estimatedCoverage) + 'x'
		if not os.path.isfile(report_file):
			writer = open(report_file, 'wt')
		else:
			writer = open(report_file, 'at')
		writer.write(report + '\n')
		writer.flush()
		writer.close()

		report = 'Estimated depth coverage: ' + str(estimatedCoverage) + 'x'
		if estimatedCoverage >= 30:
			pass_qc = True
			print report
		else:
			failing['sample'] = report + ' (lower than 30x)'
			print failing['sample']

	else:
		failing['sample'] = 'Did not run'
		print failing['sample']

	return run_successfully, pass_qc, failing
