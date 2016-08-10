import utils
import os
import multiprocessing


def fastQintegrity(fastq, outdir):
	run_successfully = False

	temporary_output_file = os.path.join(outdir, os.path.splitext(os.path.basename(fastq))[0])

	command = ['', '--stdout', '--keep', fastq, '>', temporary_output_file]

	filetype = utils.compressionType(fastq)
	if filetype == 'gz':
		command[0] = 'gunzip'
	elif filetype == 'bz2':
		command[0] = 'bunzip2'

	if command[0] != '':
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None)

	if os.path.isfile(temporary_output_file):
		os.remove(temporary_output_file)

	utils.saveVariableToPickle(run_successfully, outdir, os.path.basename(fastq))


# Count sequenced bases
def runFastQintegrity(fastq_files, threads, outdir):
	failing = {}
	failing['sample'] = False
	not_corruption_found = True

	# failing[reads].append('Bad per base N content')

	# Create Trimmomatic output directory
	fastQintegrity_folder = os.path.join(outdir, 'fastQintegrity', '')
	utils.removeDirectory(fastQintegrity_folder)
	os.mkdir(fastQintegrity_folder)

	pool = multiprocessing.Pool(processes=threads)
	for fastq in fastq_files:
		pool.apply_async(fastQintegrity, args=(fastq, fastQintegrity_folder,))
	pool.close()
	pool.join()

	files = [f for f in os.listdir(fastQintegrity_folder) if not f.startswith('.') and os.path.isfile(os.path.join(fastQintegrity_folder, f))]
	for file_found in files:
		if file_found.endswith('.pkl'):
			file_run_successfully = utils.extractVariableFromPickle(os.path.join(fastQintegrity_folder, file_found))
			if not file_run_successfully:
				failing[os.path.splitext(file_found)[0]] = ['The file is possibly corrupt']
				print os.path.splitext(file_found)[0] + ': the file is possibly corrupt'
		os.remove(os.path.join(fastQintegrity_folder, file_found))

	if len(failing) > 1:
		failing.pop('sample')
		not_corruption_found = False

	utils.removeDirectory(fastQintegrity_folder)

	return not_corruption_found, failing
