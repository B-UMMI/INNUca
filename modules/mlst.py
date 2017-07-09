import utils
import os
from functools import partial
import sys


mlst_timer = partial(utils.timer, name='MLST')


def get_species_scheme_map_version(mlst_folder):
    species_scheme_map_version = 1

    mlst_db_path = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'species_scheme_map.tab')
    if not os.path.isfile(mlst_db_path):
        mlst_db_path = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'scheme_species_map.tab')
        if not os.path.isfile(mlst_db_path):
            sys.exit('ERROR: species_scheme_map not found. Contact the developers. In the meantime try running INNUca with --skipMLST option')
        else:
            species_scheme_map_version = 2
    return mlst_db_path, species_scheme_map_version


def set_species_scheme_map_variables(list_values, species_scheme_map_version):
    if species_scheme_map_version == 1:
        val_genus = list_values[0]
        val_species = list_values[1]
        val_scheme = list_values[2]
    elif species_scheme_map_version == 2:
        val_genus = list_values[1]
        val_species = list_values[2]
        val_scheme = list_values[0]
    return val_genus, val_species, val_scheme


def parse_species_scheme_map(species_splited, mlst_db_path, species_scheme_map_version):
    scheme = 'unknown'
    genus_mlst_scheme = None

    with open(mlst_db_path, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.lower().split('\t')
                    line = [line[i].split(' ')[0] for i in range(0, len(line))]
                    val_genus, val_species, val_scheme = set_species_scheme_map_variables(line, species_scheme_map_version)
                    if val_genus == species_splited[0]:
                        if val_species == '':
                            genus_mlst_scheme = val_scheme
                        elif val_species == species_splited[1]:
                            scheme = val_scheme
        if scheme == 'unknown' and genus_mlst_scheme is not None:
            scheme = genus_mlst_scheme
    return scheme, genus_mlst_scheme


def getScheme(species):
    command = ['which', 'mlst']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    mlst_folder = os.path.abspath(os.path.realpath(stdout.splitlines()[0]))

    mlst_db_path, species_scheme_map_new = get_species_scheme_map_version(mlst_folder)

    scheme, genus_mlst_scheme = parse_species_scheme_map(species.lower().split(' '), mlst_db_path, species_scheme_map_new)

    print '\n' + 'MLST scheme found for {species}: {scheme}'.format(species=species, scheme=scheme)

    return scheme, species.lower().split(' ')[0], genus_mlst_scheme


def getBlastPath():
    print '\n' + 'The following blastn will be used'
    command = ['which', 'blastn']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    print stdout


@mlst_timer
def runMlst(contigs, scheme, outdir, species_genus, mlst_scheme_genus):
    pass_qc = False
    failing = {}
    failing['sample'] = False
    warnings = {}

    command = ['mlst', contigs]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    if run_successfully:
        scheme_mlst = stdout.splitlines()[0].split('\t')[1].split('_')[0]
        st = stdout.splitlines()[0].split('\t')[2]
        profile = stdout.splitlines()[0].split('\t')[3:]

        report = 'MLST found ST ' + str(st) + ' from scheme ' + scheme_mlst
        print report
        with open(os.path.join(outdir, 'mlst_report.txt'), 'wt') as writer:
            writer.write('#scheme' + '\n' + scheme_mlst + '\n' + '#ST' + '\n' + st + '\n')
            writer.write('#profile' + '\n' + ' '.join(profile) + '\n')
            writer.flush()

        if scheme_mlst.split('_', 1)[0] == scheme.split('_', 1)[0]:
            pass_qc = True
        else:
            if scheme == 'unknown' and scheme_mlst != '-':
                pass_qc = True
                warnings['sample'] = 'Found {scheme_mlst} scheme for a species with unknown scheme'.format(scheme_mlst=scheme_mlst)
            elif scheme == 'unknown' and scheme_mlst == '-':
                pass_qc = True
            elif species_genus == 'yersinia' and mlst_scheme_genus == 'yersinia':
                pass_qc = True
                warnings['sample'] = 'Found a Yersinia scheme ({scheme_mlst}), but it is different from what it was expected({scheme})'.format(scheme_mlst=scheme_mlst, scheme=scheme)
            else:
                failing['sample'] = 'MLST scheme found (' + scheme_mlst + ') and provided (' + scheme + ') are not the same'
                print failing['sample']
    else:
        warnings['sample'] = 'Did not run;'
        pass_qc = True

    if len(warnings) > 0:
        print warnings['sample']

    return run_successfully, pass_qc, failing, warnings
