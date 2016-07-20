import subprocess
import time
import sys
import argparse
import os
import shutil
import shlex


def parseArguments(version):
	parser = argparse.ArgumentParser(prog='INNUca.py', description='INNUca - Reads Control and Assembly', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required_options = parser.add_argument_group('Required options')
	required_options.add_argument('-i', '--inputDirectory', nargs=1, type=str, metavar='/path/to/input/directory/', help='Path to directory containing the fastq files. Can be organized in separete directories by samples or all together', required=True)
	required_options.add_argument('-s', '--speciesExpected', nargs=1, type=str, metavar='"Streptococcus agalactiae"', help='Expected species name', required=True)
	required_options.add_argument('-g', '--genomeSizeExpectedMb', nargs=1, type=float, metavar='2.1', help='Expected genome size in Mb', required=True)

	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	general_options = parser.add_argument_group('General options')
	general_options.add_argument('-o', '--outdir', nargs=1, type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')
	general_options.add_argument('-j', '--threads', nargs=1, type=int, metavar='N', help='Number of threads', required=False, default=[1])
	general_options.add_argument('--doNotUseProvidedSoftware', action='store_true', help='Tells the software to not use FastQC, Trimmomatic, SPAdes and Samtools that are provided with INNUca.py')
	general_options.add_argument('--pairEnd_filesSeparation', nargs=2, type=str, metavar='"_left/rigth.fq.gz"', help='For unusual pair-end files separation designations, you can provide two strings containning the end of fastq files names to designate each file from a pair-end data ("_left.fq.gz" "_rigth.fq.gz" for sample_left.fq.gz sample_right.fq.gz)', required=False, default=None)
	general_options.add_argument('--skipEstimatedCoverage', action='store_true', help='Tells the programme to not estimate coverage depth based on number of sequenced nucleotides and expected genome size')
	general_options.add_argument('--skipFastQC', action='store_true', help='Tells the programme to not run FastQC analysis')
	general_options.add_argument('--skipTrimmomatic', action='store_true', help='Tells the programme to not run Trimmomatic')
	general_options.add_argument('--skipSPAdes', action='store_true', help='Tells the programme to not run SPAdes and consequently MLST analysis (requires SPAdes contigs)')
	general_options.add_argument('--skipMLST', action='store_true', help='Tells the programme to not run MLST analysis')

	adapters_options = parser.add_mutually_exclusive_group()
	adapters_options.add_argument('--adapters', nargs=1, type=argparse.FileType('r'), metavar='adaptersFile.fasta', help='Fasta file containing adapters sequences to be used in FastQC and Trimmomatic', required=False, default=[None])
	adapters_options.add_argument('--doNotSearchAdapters', action='store_true', help='Tells INNUca.py to not search for adapters and clip them during Trimmomatic step')

	trimmomatic_options = parser.add_argument_group('Trimmomatic options')
	trimmomatic_options.add_argument('--doNotTrimCrops', action='store_true', help='Tells INNUca.py to not cut the beginning and end of reads during Trimmomatic step (unless specified with --trimCrop or --trimHeadCrop, INNUca.py will search for nucleotide content bias at both ends and will cut by there)')
	trimmomatic_options.add_argument('--trimCrop', nargs=1, type=int, metavar='N', help='Cut the specified number of bases to the end of the maximum reads length', required=False)
	trimmomatic_options.add_argument('--trimHeadCrop', nargs=1, type=int, metavar='N', help='Trimmomatic: cut the specified number of bases from the start of the reads', required=False)
	trimmomatic_options.add_argument('--trimSlidingWindow', nargs=1, type=str, metavar='window:meanQuality', help='Trimmomatic: perform a sliding window trimming, cutting once the average quality within the window falls below a threshold', required=False, default=['5:20'])
	trimmomatic_options.add_argument('--trimLeading', nargs=1, type=int, metavar='N', help='Trimmomatic: cut bases off the start of a read, if below a threshold quality', required=False, default=[3])
	trimmomatic_options.add_argument('--trimTrailing', nargs=1, type=int, metavar='N', help='Trimmomatic: cut bases off the end of a read, if below a threshold quality', required=False, default=[3])
	trimmomatic_options.add_argument('--trimMinLength', nargs=1, type=int, metavar='N', help='Trimmomatic: drop the read if it is below a specified length', required=False, default=[55])

	spades_options = parser.add_argument_group('SPAdes options')
	spades_options.add_argument('--spadesNotUseCareful', action='store_true', help='Tells SPAdes to only perform the assembly without the --careful option')
	spades_options.add_argument('--spadesMinContigsLength', nargs=1, type=int, metavar='N', help='Filter SPAdes contigs for length greater or equal than this value', required=False, default=[200])
	spades_options.add_argument('--spadesKmers', nargs=1, type=spades_kmers, metavar='55,77', help='Manually sets SPAdes k-mers lengths (all values must be odd, less than 128)', required=False, default=[55, 77, 99, 113, 127])
	spades_options.add_argument('--spadesMaxMemory', nargs=1, type=int, metavar='N', help='The maximum amount of RAM Gb for SPAdes to use', required=False, default=[25])
	spades_options.add_argument('--spadesMinCoverage', nargs=1, type=spades_cov_cutoff, metavar='10', help='The minimum number of reads to consider an edge in the de Bruijn graph (or path I am not sure). Can also be auto or off', required=False, default=['off'])

	args = parser.parse_args()

	if args.doNotTrimCrops and (args.trimCrop or args.trimHeadCrop):
		parser.error('Cannot use --doNotTrimCrops option with --trimCrop or --trimHeadCrop')

	return args


# For parseArguments
def spades_kmers(arguments):
	arguments = map(int, arguments.split(','))
	arguments = sorted(arguments)
	for number in arguments:
		if number % 2 != 0:
			if number < 128:
				continue
			else:
				argparse.ArgumentParser.error()
		else:
			argparse.ArgumentParser.error()
	return arguments


# For parseArguments
def spades_cov_cutoff(argument):
	string_options = ['auto', 'off']
	for option in string_options:
		if str(argument) == option:
			return str(argument)

	try:
		argument = int(argument)
		if argument > 0:
			return argument
		else:
			argparse.ArgumentParser.error('--spadesMinCoverage must be positive integer, auto or off')
	except:
		argparse.ArgumentParser.error('--spadesMinCoverage must be positive integer, auto or off')


def runCommandPopenCommunicate(command):
	run_successfully = False
	if isinstance(command, basestring):
		command = shlex.split(command)
	else:
		command = shlex.split(' '.join(command))
	print 'Running: ' + ' '.join(command)
	proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


def runCommandPopenCommunicateShell(command):
	run_successfully = False
	if not isinstance(command, basestring):
		command = ' '.join(command)
	print 'Running: ' + command
	proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = proc.communicate()
	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


def runTime(start_time):
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print 'Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's'
	return time_taken


class Logger(object):
	def __init__(self, out_directory, time_str):
		self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()

	def flush(self):
		pass

	def getLogFile(self):
		return self.logfile


# Check programs versions
def checkPrograms(programs_version_dictionary):
	print '\n' + 'Checking dependencies...'
	programs = programs_version_dictionary
	which_program = ['which', '']
	listMissings = []
	for program in programs:
		which_program[1] = program
		run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program)
		if not run_successfully:
			listMissings.append(program + ' not found in PATH.')
		else:
			if programs[program][0] is None:
				print program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0]
			else:
				check_version = [stdout.splitlines()[0], programs[program][0]]
				run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version)
				if stdout == '':
					stdout = stderr
				if program == 'bunzip2':
					version_line = stdout.splitlines()[0].rsplit(',', 1)[0].split(' ')[-1]
				else:
					version_line = stdout.splitlines()[0].split(' ')[-1]
				replace_characters = ['"', 'v', 'V', '+']
				for i in replace_characters:
					version_line = version_line.replace(i, '')
				print program + ' (' + version_line + ') found'
				if programs[program][1] == '>=':
					program_found_version = version_line.split('.')
					program_version_required = programs[program][2].split('.')
					if float('.'.join(program_found_version[0:2])) < float('.'.join(program_version_required[0:2])):
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
					elif float('.'.join(program_found_version[0:2])) == float('.'.join(program_version_required[0:2])):
						if len(program_version_required) == 3:
							if len(program_found_version) == 2:
								program_found_version.append(0)
							if program_found_version[2].split('_')[0] < program_version_required[2]:
								listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
				else:
					if version_line != programs[program][2]:
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
	return listMissings


def setPATHvariable(doNotUseProvidedSoftware, script_path):
	path_variable = os.environ['PATH']
	script_folder = os.path.dirname(script_path)
	# Set path to use provided softwares
	if not doNotUseProvidedSoftware:
		fastQC = os.path.join(script_folder, 'src', 'fastqc_v0.11.5')
		trimmomatic = os.path.join(script_folder, 'src', 'Trimmomatic-0.36')
		spades = os.path.join(script_folder, 'src', 'SPAdes-3.7.1-Linux', 'bin')

		os.environ['PATH'] = str(':'.join([fastQC, trimmomatic, spades, path_variable]))
	print '\n' + 'PATH variable:'
	print os.environ['PATH']


# Search for fastq files in the directory
# pairEnd_filesSeparation_list can be None
def searchFastqFiles(directory, pairEnd_filesSeparation_list, listAllFilesOnly):
	filesExtensions = ['fastq.gz', 'fq.gz', 'fastq.bz2', 'fq.bz2']
	pairEnd_filesSeparation = [['_R1_001.f', '_R2_001.f'], ['_1.f', '_2.f']]

	if pairEnd_filesSeparation_list is not None:
		filesExtensions = pairEnd_filesSeparation_list

	files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]

	if listAllFilesOnly:
		filesToUse = files
	else:
		filesToUse = []
		if len(files) >= 1:
			filesCorrectExtensions = []
			for extension in filesExtensions:
				for file in files:
					if file.endswith(extension):
						filesCorrectExtensions.append(file)
				if pairEnd_filesSeparation_list is None:
					if len(filesCorrectExtensions) >= 1:
						if len(filesCorrectExtensions) > 2:
							file_pair = []
							for PE_separation in pairEnd_filesSeparation:
								for file in filesCorrectExtensions:
									if PE_separation[0] in file or PE_separation[1] in file:
										file_pair.append(file)
								if len(file_pair) >= 1:
									break
							filesCorrectExtensions = file_pair
						break
			for file in filesCorrectExtensions:
				filesToUse.append(os.path.join(directory, file))
	return filesToUse


# Organize fastq files from multiple samples into samples folders
def organizeSamplesFastq(directory, pairEnd_filesSeparation_list):
	# Get files in directory
	files = searchFastqFiles(directory, None, True)

	if len(files) == 0:
		sys.exit('No fastq files were found!')

	pairEnd_filesSeparation = [['_R1_001.f', '_R2_001.f'], ['_1.f', '_2.f']]
	if pairEnd_filesSeparation_list is not None:
		pairEnd_filesSeparation = pairEnd_filesSeparation_list

	# Store samples files information
	samples = {}
	for fastq in files:
		sample = None
		for PE_separation in pairEnd_filesSeparation:
			if fastq.find(PE_separation[0]) > -1:
				sample = fastq[:fastq.find(PE_separation[0])]
				break
			if fastq.find(PE_separation[1]) > -1:
				sample = fastq[:fastq.find(PE_separation[1])]
				break
		if sample is not None and sample not in samples:
			samples[sample] = []
		if sample in samples:
			samples[sample].append(fastq)

	# Create the file structure required
	for sample in samples:
		if len(samples[sample]) == 1:
			print 'Only one fastq file was found: ' + str(samples[sample])
			print 'Pair-End sequencing is required. This sample will be ignored'
			samples.remove(sample)
			continue
		sample_folder = os.path.join(directory, sample, '')
		if not os.path.isdir(sample_folder):
			os.makedirs(sample_folder)
		for fastq in samples[sample]:
			link_path = os.path.join(sample_folder, os.path.basename(fastq))
			if os.path.islink(link_path):
				os.remove(link_path)
			if not os.path.isfile(link_path):
				os.symlink(os.path.join(directory, fastq), link_path)

	return samples.keys()


# Check if input directory exists with fastq files and store samples name that have fastq files
def checkSetInputDirectory(inputDirectory, outdir, pairEnd_filesSeparation_list):
	samples = []
	removeCreatedSamplesDirectories = False

	if not os.path.isdir(inputDirectory):
		sys.exit('Input directory does not exist!')
	else:
		directories = [d for d in os.listdir(inputDirectory) if not d.startswith('.') and os.path.isdir(os.path.join(inputDirectory, d, ''))]
		if os.path.basename(outdir) in directories:
			print 'Output directory is inside input directory and will be ignore in the checking and setting input directory step'
			directories.remove(os.path.basename(outdir))
		if len(directories) == 0:
			print('There is no samples folders! Search for fastq files in input directory')
			samples = organizeSamplesFastq(inputDirectory, pairEnd_filesSeparation_list)
			removeCreatedSamplesDirectories = True
		else:
			for directory in directories:
				files = searchFastqFiles(os.path.join(inputDirectory, directory, ''), pairEnd_filesSeparation_list, False)
				if len(files) == 1:
					print 'Only one fastq file was found: ' + str(files)
					print 'Pair-End sequencing is required. This sample will be ignored'
					continue
				elif len(files) >= 1:
					samples.append(directory)
	if len(samples) == 0:
		sys.exit('There is no fastq files for the samples folders provided!')
	return samples, removeCreatedSamplesDirectories


# Remove directory
def removeDirectory(directory):
	if os.path.isdir(directory):
		shutil.rmtree(directory)


# Get script version
def scriptVersionGit(version, directory, script_path):
	print 'Version ' + version
	os.chdir(os.path.dirname(script_path))
	command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
	proc = subprocess.Popen(shlex.split(' '.join(command)), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	print stdout
	command = ['git', 'remote', 'show', 'origin']
	proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	print stdout
	os.chdir(directory)
