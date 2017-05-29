import utils
import os
from functools import partial

mlst_timer = partial(utils.timer, name='MLST')


def getScheme(species):
	command = ['which', 'mlst']
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
	mlst_folder = os.path.abspath(os.path.realpath(stdout.splitlines()[0]))
	mlst_db_path = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'species_scheme_map.tab')

	species = species.lower().split(' ')

	scheme = 'unknown'

	with open(mlst_db_path, 'rtU') as reader:
		genus_mlst_scheme = None
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not line.startswith('#'):
					line = line.lower().split('\t')
					line = [line[i].split(' ')[0] for i in range(0, len(line))]
					if line[0] == species[0]:
						if line[1] == '':
							genus_mlst_scheme = line[2]
						elif line[1] == species[1]:
							scheme = line[2]
		if scheme == 'unknown' and genus_mlst_scheme is not None:
			scheme = genus_mlst_scheme

	print '\n' + 'MLST scheme found for ' + ' '.join(species) + ': ' + scheme

	return scheme


def getBlastPath():
	print '\n' + 'The following blastn will be used'
	command = ['which', 'blastn']
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	print stdout


@mlst_timer
def runMlst(contigs, scheme, outdir):
	pass_qc = False
	failing = {}
	failing['sample'] = False

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

		if scheme_mlst == scheme:
			pass_qc = True
		else:
			failing['sample'] = 'MLST scheme found (' + scheme_mlst + ') and provided (' + scheme + ') are not the same'
			print failing['sample']
	else:
		failing['sample'] = 'Did not run;'
		print failing['sample']

	return run_successfully, pass_qc, failing
