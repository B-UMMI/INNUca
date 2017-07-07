import os.path
import utils
import functools
import sys


def check_existing_default_config(species, script_path):
    species = species.lower().split(' ')
    trueCoverage_config_folder = os.path.join(os.path.dirname(script_path), 'modules', 'trueCoverage_rematch', '')
    config = None
    reference = None
    files = [f for f in os.listdir(trueCoverage_config_folder) if not f.startswith('.') and os.path.isfile(os.path.join(trueCoverage_config_folder, f))]
    for file_found in files:
        file_path = os.path.join(trueCoverage_config_folder, file_found)
        if file_found == '_'.join(species) + '.config':
            config = file_path
        elif file_found == '_'.join(species) + '.fasta':
            reference = file_path
    return config, reference


def parse_config(config_file):
    config = {'reference_file': None, 'length_extra_seq': None, 'maximum_number_absent_genes': None, 'maximum_number_genes_multiple_alleles': None, 'minimum_read_coverage': None, 'minimum_depth_presence': None, 'minimum_depth_call': None, 'minimum_depth_frequency_dominant_allele': None, 'minimum_gene_coverage': None, 'minimum_gene_identity': None}

    with open(config_file, 'rtU') as reader:
        field = None
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                line = line.split(' ')[0]
                if line.startswith('#'):
                    line = line[1:].split(' ')[0]
                    field = line
                else:
                    if field is not None:
                        if field in ['length_extra_seq', 'maximum_number_absent_genes', 'maximum_number_genes_multiple_alleles', 'minimum_read_coverage', 'minimum_depth_presence', 'minimum_depth_call', 'minimum_gene_coverage', 'minimum_gene_identity']:
                            line = int(line)
                            if field in ['minimum_gene_coverage', 'minimum_gene_identity']:
                                if line < 0 or line > 100:
                                    sys.exit('minimum_gene_coverage in trueCoverage_rematch config file must be an integer between 0 and 100')
                        elif field == 'minimum_depth_frequency_dominant_allele':
                            line = float(line)
                            if line < 0 or line > 1:
                                sys.exit('minimum_depth_frequency_dominant_allele in trueCoverage_rematch config file must be a double between 0 and 1')
                        config[field] = line
                        field = None

    for field in config:
        if config[field] is None:
            sys.exit(field + ' in trueCoverage_rematch config file is missing')

    return config


def clean_headers_reference_file(reference_file, outdir, extraSeq, rematch_module):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    print 'Checking if reference sequences contain ' + str(problematic_characters) + '\n'
    headers_changed = False
    new_reference_file = str(reference_file)
    sequences, genes, headers_changed = rematch_module.get_sequence_information(reference_file, extraSeq)
    if headers_changed:
        print 'At least one of the those characters was found. Replacing those with _' + '\n'
        new_reference_file = os.path.join(outdir, os.path.splitext(os.path.basename(reference_file))[0] + '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = rematch_module.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, genes, sequences


trueCoverage_timer = functools.partial(utils.timer, name='True coverage check')


@trueCoverage_timer
def runTrueCoverage(sample, fastq, reference, threads, outdir, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, conserved_True, debug, numMapLoc, minGeneIdentity, trueCoverage_config, rematch_script):
    pass_qc = False
    failing = {}

    trueCoverage_folder = os.path.join(outdir, 'trueCoverage', '')
    utils.removeDirectory(trueCoverage_folder)
    os.mkdir(trueCoverage_folder)

    sys.path.append(os.path.join(os.path.dirname(rematch_script), 'modules', ''))
    import rematch_module

    # Run ReMatCh
    reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(reference, trueCoverage_folder, extraSeq, rematch_module)
    time_taken, run_successfully, data_by_gene, sample_data_general, consensus_files, consensus_sequences = rematch_module.runRematchModule(sample, fastq, reference_file, threads, trueCoverage_folder, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, True, debug, 1, minGeneIdentity, 'first', 7, 'none', reference_dict, 'X', None, gene_list_reference, True)

    if run_successfully:
        print 'Writing report file'
        os.rename(os.path.join(trueCoverage_folder, 'rematchModule_report.txt'), os.path.join(outdir, 'trueCoverage_report.txt'))

        if sample_data_general['number_absent_genes'] > trueCoverage_config['maximum_number_absent_genes']:
            failing['absent_genes'] = 'The number of absent genes (' + str(sample_data_general['number_absent_genes']) + ') exceeds the maximum allowed (' + str(trueCoverage_config['maximum_number_absent_genes']) + ')'
        if sample_data_general['number_genes_multiple_alleles'] > trueCoverage_config['maximum_number_genes_multiple_alleles']:
            failing['multiple_alleles'] = 'The number of genes with multiple alleles (' + str(sample_data_general['number_genes_multiple_alleles']) + ') exceeds the maximum allowed (' + str(trueCoverage_config['maximum_number_genes_multiple_alleles']) + ')'
        if sample_data_general['mean_sample_coverage'] < trueCoverage_config['minimum_read_coverage']:
            failing['read_coverage'] = 'The mean read coverage for genes present (' + str(sample_data_general['mean_sample_coverage']) + ') dit not meet the minimum required (' + str(trueCoverage_config['minimum_read_coverage']) + ')'
    else:
        failing['sample'] = 'Did not run'

    if len(failing) == 0:
        pass_qc = True
        failing['sample'] = False
    else:
        print failing

    if not debug:
        utils.removeDirectory(trueCoverage_folder)

    return run_successfully, pass_qc, failing
