import utils
import os
from functools import partial

mlst_timer = partial(utils.timer, name='MLST')


@mlst_timer
def runMlst(contigs, species, outdir):
	pass_qc = False
	failing = {}
	failing['sample'] = False

	species = species.lower().split(' ')
	species = species[0][0] + species[1]

	command = ['mlst', contigs]

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	if run_successfully:
		species_mlst = stdout.splitlines()[0].split('\t')[1].split('_')[0]
		st = stdout.splitlines()[0].split('\t')[2]

		report = 'MLST found ST ' + str(st) + ' from species ' + species_mlst
		print report
		with open(os.path.join(outdir, 'mlst_report.txt'), 'wt') as writer:
			writer.write('#species' + '\n' + species_mlst + '\n' + '#ST' + '\n' + st + '\n')
			writer.flush()

		if species_mlst == species:
			pass_qc = True
		else:
			failing['sample'] = 'Species found (' + species_mlst + ') and provided (' + species + ') are not the same'
			print failing['sample']
	else:
		failing['sample'] = 'Did not run;'
		print failing['sample']

	return run_successfully, pass_qc, failing
