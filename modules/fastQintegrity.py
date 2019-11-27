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

    encoding, min_reads_length, max_reads_length, num_reads, num_bp = None, None, None, None, None

    if compression_type is not None:
        command = [compression_type[1], '--stdout', '--keep', fastq, '>', temporary_output_file]
        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)

        if run_successfully:
            encoding, min_reads_length, max_reads_length, num_reads, num_bp = \
                run_guess_encoding_single_thread(temporary_output_file, None, outdir)

    if os.path.isfile(temporary_output_file):
        os.remove(temporary_output_file)

    utils.saveVariableToPickle([run_successfully, encoding, min_reads_length, max_reads_length, num_reads, num_bp],
                               outdir, os.path.basename(fastq))


def run_guess_encoding_single_thread(fastq_file, number_reads_access_None_all, outdir):
    outdir_guess_encoding = os.path.join(outdir, os.path.splitext(os.path.basename(fastq_file))[0])
    utils.removeDirectory(outdir_guess_encoding)
    os.mkdir(outdir_guess_encoding)

    guess_encoding.guess_encoding(fastq_file, number_reads_access_None_all, outdir_guess_encoding)
    encoding_data = guess_encoding.gather_data_together(outdir_guess_encoding)
    final_enconding = guess_encoding.get_final_encoding(encoding_data)

    min_reads_length, max_reads_length, _, _ = guess_encoding.determine_min_max_reads_length(encoding_data)
    num_reads, num_bp = guess_encoding.get_num_reads_bp(encoding_data)

    utils.removeDirectory(outdir_guess_encoding)
    return final_enconding, min_reads_length, max_reads_length, num_reads, num_bp


def report_reads_length(min_reads_length_each_fastq, max_reads_length_each_fastq, outdir):
    """
    Writes reads length report

    Parameters
    ----------
    min_reads_length_each_fastq : list
        Minimum reads length found for each fastq file
    max_reads_length_each_fastq : list
        Maximum reads length found for each fastq file
    outdir : str
        Path to the output directory

    Returns
    -------

    """

    with open(os.path.join(outdir, 'reads_length_report.tab'), 'wt') as writer:
        writer.write('#' + '\t'.join(['min', 'max']) + '\n')
        writer.write('\t'.join([';'.join(map(str, set(min_reads_length_each_fastq))),
                                ';'.join(map(str, set(max_reads_length_each_fastq)))]) + '\n')


def report_num_reads_bp(num_reads, num_bp, outdir):
    """
    Writes the total number of reads and bp sequenced

    Parameters
    ----------
    num_reads : int
        Total number of reads sequenced
    num_bp : int
        Total number of bp sequenced
    outdir : str
        Path to the output directory

    Returns
    -------

    """

    with open(os.path.join(outdir, 'num_reads_bp_report.tab'), 'wt') as writer:
        writer.write('#' + '\t'.join(['num_reads', 'num_bp']) + '\n')
        writer.write('\t'.join(map(str, [num_reads, num_bp])) + '\n')


fastq_timer = partial(utils.timer, name='FastQ integrity check')


@fastq_timer
def runFastQintegrity(fastq_files, threads, outdir):
    pass_qc = True
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

    encoding = {}
    num_reads, num_bp = 0, 0
    files = [f for f in os.listdir(fastQintegrity_folder) if not f.startswith('.') and os.path.isfile(os.path.join(fastQintegrity_folder, f))]
    for file_found in files:
        if file_found.endswith('.pkl'):
            file_run_successfully, file_encoding, min_reads_length, max_reads_length, num_reads_fastq, num_bp_fastq = \
                utils.extractVariableFromPickle(os.path.join(fastQintegrity_folder, file_found))
            if file_run_successfully:
                encoding[file_found] = {'file_encoding': file_encoding,
                                        'min_reads_length': min_reads_length,
                                        'max_reads_length': max_reads_length}
                num_reads += num_reads_fastq if num_reads_fastq is not None else 0
                num_bp += num_bp_fastq if num_bp_fastq is not None else 0
            else:
                failing[os.path.splitext(file_found)[0]] = ['The file is possibly corrupt']
                print(os.path.splitext(file_found)[0] + ': the file is possibly corrupt')
        os.remove(os.path.join(fastQintegrity_folder, file_found))

    if len(failing) > 1:
        failing.pop('sample')
        not_corruption_found = False
        pass_qc = False

        min_reads_length_found, max_reads_length_found, num_reads, num_bp = None, None, None, None

    if len(encoding) == 0:
        encoding = None
        print('It was no possible to determine the FASTQ encodings')
    else:
        min_reads_length_found, max_reads_length_found, min_reads_length_each_fastq, max_reads_length_each_fastq = \
            guess_encoding.determine_min_max_reads_length(encoding)
        report_reads_length(min_reads_length_each_fastq, max_reads_length_each_fastq, outdir)

        if len(set([x['file_encoding'][0] for x in encoding.values() if x['file_encoding'] is not None])) == 1:
            encoding = [x['file_encoding'][0] for x in encoding.values() if x['file_encoding'] is not None][0]
            print('Fastq quality encoding: {0}'.format(str(encoding)))
        else:
            print('It was no possible to determine the FASTQ encodings')
            print('This was what has been found: {0}'.format(str(encoding)))
            encoding = None

    utils.removeDirectory(fastQintegrity_folder)

    return not_corruption_found, pass_qc, failing, encoding, min_reads_length_found, max_reads_length_found
