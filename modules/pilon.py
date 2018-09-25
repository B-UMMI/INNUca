import utils
import os
from functools import partial


# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
    if os.path.isfile(str(referenceFile + '.1.bt2')):
        run_successfully = True
    else:
        command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    return run_successfully


# Mapping with Bowtie2
def mapping_bowtie2(fastq_files, reference_file, outdir, keep_bam=False, threads=1):
    """
    Map reads against a reference fasta file

    Parameters
    ----------
    fastq_files : list
        List of fastq files
    reference_file : str
        Path to the reference file (the assembly)
    outdir : str
        Path to the output directory
    keep_bam : bool, default False
        True if want to keep the BAM file produced (with mapped and unmapped reads)
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Assembly_Mapping module ran successfully or not
    sam_file : str or None
        If everything went fine, it returns the path to the sam file, otherwise it returns None
    """

    sam_file = os.path.join(outdir, str('alignment.sam'))

    # Index reference file
    run_successfully = indexSequenceBowtie2(reference_file, threads)

    if run_successfully:
        command = ['bowtie2', '-q', '--very-sensitive-local', '--threads', str(threads), '-x', reference_file, '',
                   '', '-S', sam_file]
        if len(fastq_files) == 1:
            command[7] = '-U ' + fastq_files[0]
        else:
            command[7] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]

        if not keep_bam:
            command[8] = '--no-unal'

        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    if not run_successfully:
        sam_file = None

    return run_successfully, sam_file


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    if not run_successfully:
        output_file = None
    return run_successfully, output_file


# Index alignment file
def indexAlignment(alignment_file):
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    return run_successfully


def pilon(jar_path_pilon, assembly, bam_file, outdir, jarMaxMemory):
    assembly_polished = os.path.splitext(assembly)[0] + '.polished.fasta'
    command = ['java', '', '-jar', jar_path_pilon, '--genome', assembly, '--frags', bam_file, '--outdir', outdir, '--output', os.path.basename(os.path.splitext(assembly_polished)[0]), '--changes', '--vcf']
    if str(jarMaxMemory) != 'off':
        command[1] = '-Xmx' + str(int(round(jarMaxMemory * 1024, 0))) + 'M'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
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


pilon_timer = partial(utils.timer, name='Pilon')


@pilon_timer
def run_pilon(jar_path_pilon, assembly, fastq_files, outdir, jar_max_memory, alignment_file, keep_bam=False, threads=1):
    """
    Runs Assembly_Mapping for INNUca and QA/QC the results

    Parameters
    ----------
    jar_path_pilon
    assembly : str
        Path to the assembly to correct
    fastq_files : list
        List of fastq files
    outdir : str
        Path to the output directory
    jar_max_memory : int or 'off'
        If not 'off' is provided, sets the maximum RAM Gb usage by jar files
    alignment_file : str or None
        Path to the BAM file to be used. If None is provided, new alignment reads will be performed
    keep_bam : bool, default False
        True if want to keep the BAM file produced (with mapped and unmapped reads)
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Assembly_Mapping module ran successfully or not
    pass_qc : None
        QA/QC not performed
    time_taken : float
        Seconds that run_assembly_mapping took to run
    failing : dict
        Dictionary with the failing reasons. If sample did not fail, it is only {'sample': False}. If it failed, keys
        will be the level of failing, and values list of strings
    assembly_polished : str or None
        Path to the polished assembly. If something went wrong, None is returned
    pilon_folder : str
        Path to Pilon working directory
    new_bam : bool
        True if new alignment reads was performed
    alignment_file : str or None
        Path to the BAM file used to correct the assembly. If something went wrong, None is returned.
    """

    failing = {'sample': False}

    pilon_folder = os.path.join(outdir, 'pilon', '')
    utils.removeDirectory(pilon_folder)
    os.mkdir(pilon_folder)

    # Create a symbolic link to the assembly
    assembly_link = os.path.join(pilon_folder, os.path.basename(assembly))
    os.symlink(assembly, assembly_link)

    run_successfully = True

    new_bam = False
    if alignment_file is None:
        # Index assembly using Bowtie2
        run_successfully = indexSequenceBowtie2(assembly_link, threads)

        if run_successfully:
            # mapping_bowtie2(fastq_files, reference_file, outdir, keep_bam=False, threads=1
            run_successfully, sam_file = mapping_bowtie2(fastq_files=fastq_files, reference_file=assembly_link,
                                                         outdir=pilon_folder, keep_bam=keep_bam, threads=threads)

            if run_successfully:
                alignment_file = os.path.splitext(sam_file)[0] + '.bam'
                run_successfully, alignment_file = sortAlignment(sam_file, alignment_file, False, threads)

                if run_successfully:
                    os.remove(sam_file)
                    run_successfully = indexAlignment(alignment_file)
                    new_bam = True
                else:
                    alignment_file = None

    assembly_polished = None

    if run_successfully:
        run_successfully, assembly_polished = pilon(jar_path_pilon, assembly_link, alignment_file, pilon_folder,
                                                    jar_max_memory)

        if run_successfully:
            parsePilonResult(assembly_polished, outdir)
            os.rename(assembly_polished, os.path.join(outdir, os.path.basename(assembly_polished)))
            assembly_polished = os.path.join(outdir, os.path.basename(assembly_polished))
            if keep_bam and new_bam:
                os.rename(alignment_file, os.path.join(outdir, '{}.bam'.format(os.path.basename(assembly))))
                alignment_file = os.path.join(outdir, '{}.bam'.format(os.path.basename(assembly)))

    if alignment_file is not None and os.path.isfile(str(alignment_file)) and not keep_bam:
        os.remove(alignment_file)

    if not run_successfully:
        failing['sample'] = 'Did not run'
        print failing['sample']

    return run_successfully, None, failing, assembly_polished, pilon_folder, new_bam, alignment_file
