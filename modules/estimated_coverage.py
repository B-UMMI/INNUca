import utils
import os
from functools import partial
import multiprocessing


# Count sequenced bases
def countSequencedBases(fastq_file, outdir):
	run_successfully = False
	bases = None

	command = ['', '--keep', '--stdout', fastq_file, '|', 'grep', '--after-context=1', '"@"', '|', 'grep', '--invert-match', '"^--$"', '|', 'grep', '--invert-match', '"@"', '|', 'wc', '']

	# Determine compression type
	filetype = utils.compressionType(fastq_file)
	if filetype != 'gz' and filetype != 'bz2':
		utils.saveVariableToPickle([run_successfully, bases], outdir, str('estimate_coverage.' + os.path.basename(fastq_file)))
	else:
		if filetype == 'gz':
			command[0] = 'gunzip'
		else:
			command[0] = 'bunzip2'

		# Number of characters
		command[18] = '--chars'
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)
		if run_successfully:
			bases = int(stdout.splitlines()[0])

			# Number of lines
			command[18] = '--lines'
			run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)
			if run_successfully:
				lines = int(stdout.splitlines()[0])
				bases = bases - lines

		utils.saveVariableToPickle([run_successfully, bases], outdir, str('estimate_coverage.' + os.path.basename(fastq_file)))


coverage_timer = partial(utils.timer, name='Estimated Coverage analysis')


# Get estimated coverage
@coverage_timer
def getEstimatedCoverage(fastq_files, estimatedGenomeSizeMb, outdir, threads):
	run_successfully = False
	pass_qc = False
	failing = {}
	failing['sample'] = False

	# Run Estimated Coverage
	estimatedCoverage = None

	# Get number bases for each fastq file
	pool = multiprocessing.Pool(processes=threads)
	for fastq in fastq_files:
		pool.apply_async(countSequencedBases, args=(fastq, outdir,))
	pool.close()
	pool.join()

	numberBases = 0
	file_problems = False
	files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
	for file_found in files:
		if file_found.startswith('estimate_coverage.') and file_found.endswith('.pkl'):
			file_path = os.path.join(outdir, file_found)

			if not file_problems:
				run_successfully, bases = utils.extractVariableFromPickle(file_path)
				if run_successfully:
					numberBases += bases
				else:
					file_problems = True

			os.remove(file_path)

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

	return run_successfully, pass_qc, failing, estimatedCoverage
