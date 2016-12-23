import utils
import os
import multiprocessing
from functools import partial
import guess_encoding


@utils.trace_unhandled_exceptions
def fastQintegrity(fastq, outdir):
	run_successfully = False

	temporary_output_file = os.path.join(outdir, os.path.splitext(os.path.basename(fastq))[0])

	compression_type = utils.compressionType(fastq)

	encoding = None

	if compression_type is not None:
		command = [compression_type[1], '--stdout', '--keep', fastq, '>', temporary_output_file]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)

		if run_successfully:
			encoding = run_guess_encoding_single_thread(temporary_output_file, None, outdir)

	if os.path.isfile(temporary_output_file):
		os.remove(temporary_output_file)

	utils.saveVariableToPickle([run_successfully, encoding], outdir, os.path.basename(fastq))


def run_guess_encoding_single_thread(fastq_file, number_reads_access_None_all, outdir):
	guess_encoding.gess_encoding(fastq_file, number_reads_access_None_all, outdir)
	encoding_data = guess_encoding.gather_data_together(outdir)
	final_enconding = guess_encoding.get_final_encoding(encoding_data)
	return final_enconding


fastq_timer = partial(utils.timer, name='FastQ integrity check')


@fastq_timer
def runFastQintegrity(fastq_files, threads, outdir):
	failing = {}
	failing['sample'] = False
	not_corruption_found = True

	fastQintegrity_folder = os.path.join(outdir, 'fastQintegrity', '')
	utils.removeDirectory(fastQintegrity_folder)
	os.mkdir(fastQintegrity_folder)

	pool = multiprocessing.Pool(processes=threads)
	for fastq in fastq_files:
		pool.apply_async(fastQintegrity, args=(fastq, fastQintegrity_folder,))
	pool.close()
	pool.join()

	encoding = []
	files = [f for f in os.listdir(fastQintegrity_folder) if not f.startswith('.') and os.path.isfile(os.path.join(fastQintegrity_folder, f))]
	for file_found in files:
		if file_found.endswith('.pkl'):
			file_run_successfully, file_encoding = utils.extractVariableFromPickle(os.path.join(fastQintegrity_folder, file_found))
			if file_run_successfully:
				encoding.append(file_encoding)
			else:
				failing[os.path.splitext(file_found)[0]] = ['The file is possibly corrupt']
				print os.path.splitext(file_found)[0] + ': the file is possibly corrupt'
		os.remove(os.path.join(fastQintegrity_folder, file_found))

	if len(failing) > 1:
		failing.pop('sample')
		not_corruption_found = False

	if len(encoding) == 0:
		encoding = None
		print 'It was no possible to determine the FASTQ encodings'
	else:
		if len(set([x[0] for x in encoding])) == 1:
			encoding = encoding[0]
			print 'Fastq quality encoding: ' + str(encoding)
		else:
			print 'It was no possible to determine the FASTQ encodings'
			print 'This was what has been found: ' + str(encoding)
			encoding = None

	utils.removeDirectory(fastQintegrity_folder)

	return not_corruption_found, None, failing, encoding  # None added for consistency with other steps
