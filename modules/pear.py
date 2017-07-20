import os.path
import utils
import multiprocessing
import functools


def check_uncompression_fastq(file_to_test):
    read_fields = [False, False, False, False]
    length_sequence = 0
    malformated_fastq = False

    with open(file_to_test, 'rtU') as reader:
        counter = 0
        while counter < 4:
            line = reader.next().splitlines()[0]
            if len(line) > 0:
                if counter == 0 and line.startswith('@'):
                    read_fields[0] = True
                elif counter == 1:
                    line = list(line.lower())
                    length_sequence = len(line)
                    if len(set(line).difference(list('ATGCN'.lower()))) == 0:
                        read_fields[1] = True
                elif counter == 2 and line.startswith('+'):
                    read_fields[2] = True
                elif counter == 3:
                    if len(line.splitlines()[0]) == length_sequence:
                        read_fields[3] = True
            counter += 1

    if not all(read_fields):
        malformated_fastq = True

    return malformated_fastq, length_sequence


@utils.trace_unhandled_exceptions
def compress_decompress(compressed_file, decompressed_file, compressed_True):
    run_successfully = False
    malformated_fastq = False
    length_sequence = None

    compression_type = None
    if not compressed_True:
        compression_type = utils.compressionType(compressed_file)

    if compression_type is not None or compressed_True:
        command = ['', '', '--stdout', '--keep', '', '>', '']

        if not compressed_True:
            command[0] = compression_type[0]
            command[1] = '--decompress'
            command[4] = compressed_file
            command[6] = decompressed_file
        else:
            command[0] = 'gzip'
            command[4] = decompressed_file
            command[6] = compressed_file

        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)
        if run_successfully and not compressed_True:
            malformated_fastq, length_sequence = check_uncompression_fastq(decompressed_file)
    elif compression_type is None and not compressed_True:
        run_successfully = True
        malformated_fastq, length_sequence = check_uncompression_fastq(compressed_file)
        decompressed_file = compressed_file

    if malformated_fastq:
        run_successfully = False

    utils.saveVariableToPickle([run_successfully, compressed_file if compressed_True else decompressed_file, length_sequence], os.path.dirname(decompressed_file), os.path.splitext(os.path.basename(decompressed_file))[0])


def get_compressed_decompressed_reads(outdir):
    run_successfully = True
    reads = {}

    counter = 0
    files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
    for file_found in files:
        if file_found.endswith('.pkl'):
            file_path = os.path.join(outdir, file_found)

            if run_successfully:
                run_successfully, decompressed_file, length_sequence = utils.extractVariableFromPickle(file_path)
                if run_successfully:
                    reads[counter] = [decompressed_file, length_sequence]
                    counter += 1

            os.remove(file_path)

    return run_successfully, [reads[i][0] for i in reads]


def parse_pearOutput_getAssembled(stdout):
    assembled_reads = None
    unassembled_reads = None
    discarded_reads = None
    for line in stdout.splitlines():
        if len(line) > 0:
            if line.startswith('Assembled reads') and line.endswith('%)'):
                assembled_reads = int(line.split('/', 1)[0].split(':', 1)[1].strip().replace(',', ''))
            elif line.startswith('Not assembled reads') and line.endswith('%)'):
                unassembled_reads = int(line.split('/', 1)[0].split(':', 1)[1].strip().replace(',', ''))
            elif line.startswith('Discarded reads') and line.endswith('%)'):
                discarded_reads = int(line.split('/', 1)[0].split(':', 1)[1].strip().replace(',', ''))

    return assembled_reads, unassembled_reads, discarded_reads


def run_pear(decompressed_reads_list, sample_name, threads, outdir, fastq_encoding, trimmomatic_run_successfully, minimum_overlap_reads):
    pass_qc = False
    failing = {}

    command = ['pear', '--forward-fastq', decompressed_reads_list[0], '--reverse-fastq', decompressed_reads_list[1], '--output', os.path.join(outdir, sample_name), '--p-value', str(1.0), '--min-assembly-length', str(minimum_overlap_reads), '--phred-base', '', '--cap', str(0), '--threads', str(threads), '--memory', str(str(threads) + 'G'), '--keep-original']

    if trimmomatic_run_successfully:
        command[12] = '33'
        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)
    else:
        if fastq_encoding is not None:
            command[12] = str(fastq_encoding[1])
            run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
        else:
            print 'Pear fail! Trying run with Phred+33 enconding defined...'
            command[12] = '33'
            run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
            if not run_successfully:
                print 'Pear fail again! Trying run with Phred+64 enconding defined...'
                command[12] = '64'
                run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    with open(os.path.join(outdir, str(sample_name + '.pear_out.txt')), 'wt') as writer:
        for line in stdout:
            writer.write(line)

    unassembled_pe_reads, assembled_se_reads = None, None
    assembled_reads, unassembled_reads, discarded_reads = None, None, None
    if run_successfully:
        assembled_reads, unassembled_reads, discarded_reads = parse_pearOutput_getAssembled(stdout)
        unassembled_pe_reads_uncompressed, assembled_se_reads_uncompressed = get_pear_output(outdir, sample_name)
        if assembled_reads == 0:
            assembled_se_reads = None
        else:
            compress_decompress(str(assembled_se_reads_uncompressed + '.gz'), assembled_se_reads_uncompressed, True)
            run_successfully, assembled_se_reads = get_compressed_decompressed_reads(outdir)
            assembled_se_reads = assembled_se_reads[0] if run_successfully else assembled_se_reads

        if unassembled_reads == 0:
            unassembled_pe_reads = None
        else:
            if run_successfully:
                pool = multiprocessing.Pool(processes=threads)
                for fastq in unassembled_pe_reads_uncompressed:
                    pool.apply_async(compress_decompress, args=(str(fastq + '.gz'), fastq, True,))
                pool.close()
                pool.join()

                run_successfully, unassembled_pe_reads = get_compressed_decompressed_reads(outdir)

        os.remove(assembled_se_reads_uncompressed)
        for fastq in unassembled_pe_reads_uncompressed:
            os.remove(fastq)

        if float(assembled_reads) / float(assembled_reads + unassembled_reads) < 0.75:
            pass_qc = True
            failing['sample'] = False
        else:
            failing['sample'] = 'Number of overlapping reads is >= 75% of total reads'
            print failing

    return run_successfully, pass_qc, failing, assembled_se_reads, unassembled_pe_reads, assembled_reads, unassembled_reads, discarded_reads


def get_pear_output(outdir, sample_name):
    unassembled_pe_reads = []
    assembled_se_reads = None
    files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
    for file_found in files:
        file_path = os.path.join(outdir, file_found)
        if file_found.startswith(sample_name + '.unassembled.'):
            unassembled_pe_reads.append(file_path)
        elif file_found.startswith(sample_name + '.assembled.'):
            assembled_se_reads = file_path
        elif file_found.startswith(sample_name + '.discarded.'):
            os.remove(file_path)

    return unassembled_pe_reads, assembled_se_reads


def determine_minimum_overlap(pear_min_overlap, min_reads_legth, max_reads_legth):
    if pear_min_overlap is None:
        if max_reads_legth is not None:
            pear_min_overlap = int(round(float(max_reads_legth) / 4 * 3, 0))
            print 'The minimum reads overlaping for Pear assembly them into only one read is ' + str(pear_min_overlap) + ' nts (3/4 of ' + str(max_reads_legth) + ' maximum reads length)'
        else:
            pear_min_overlap = int(round(float(min_reads_legth) / 4 * 3, 0))
            print 'The minimum reads overlaping for Pear assembly them into only one read is ' + str(pear_min_overlap) + ' nts (3/4 of ' + str(min_reads_legth) + ' minimum reads length)'
    return pear_min_overlap


pear_timer = functools.partial(utils.timer, name='Pear')


@pear_timer
def runPear(fastq_files, threads, outdir, sampleName, fastq_encoding, trimmomatic_run_successfully, minimum_overlap_reads):
    failing = {'sample': False}
    warnings = {}

    pear_folder = os.path.join(outdir, 'pear', '')
    utils.removeDirectory(pear_folder)
    os.mkdir(pear_folder)

    pool = multiprocessing.Pool(processes=threads)
    for fastq in fastq_files:
        pool.apply_async(compress_decompress, args=(fastq, os.path.join(pear_folder, str('temp.' + os.path.splitext(os.path.basename(fastq))[0])), False,))
    pool.close()
    pool.join()

    run_successfully, decompressed_reads = get_compressed_decompressed_reads(pear_folder)

    assembled_se_reads = None
    unassembled_pe_reads = None
    if run_successfully:
        if len(decompressed_reads) == 2:
            run_successfully, pass_qc, warnings, assembled_se_reads, unassembled_pe_reads, assembled_reads, unassembled_reads, discarded_reads = run_pear(decompressed_reads, sampleName, threads, pear_folder, fastq_encoding, trimmomatic_run_successfully, minimum_overlap_reads)
            if warnings['sample'] is False:
                warnings = {}
            if run_successfully:
                with open(os.path.join(outdir, str('pear_report.txt')), 'wt') as writer:
                    writer.write('#assembled_reads' + '\n' + str(assembled_reads) + '\n')
                    writer.write('#unassembled_reads' + '\n' + str(unassembled_reads) + '\n')
                    writer.write('#discarded_reads' + '\n' + str(discarded_reads) + '\n')
        else:
            run_successfully = False

        for fastq in decompressed_reads:
            os.remove(fastq)

    if not run_successfully:
        warnings['sample'] = 'Did not run'
        print warnings

    return run_successfully, True, failing, unassembled_pe_reads, assembled_se_reads, pear_folder, warnings
