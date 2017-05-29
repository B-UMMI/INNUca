import subprocess
import time
import sys
import argparse
import os
import shutil
import shlex
from threading import Timer
import pickle
import functools
import traceback
import csv


def parseArguments(version):
	parser = argparse.ArgumentParser(prog='INNUca.py', description='INNUca - Reads Control and Assembly', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required_options = parser.add_argument_group('Required options')
	required_options.add_argument('-s', '--speciesExpected', type=str, metavar='"Streptococcus agalactiae"', help='Expected species name', required=True)
	required_options.add_argument('-g', '--genomeSizeExpectedMb', type=float, metavar='2.1', help='Expected genome size in Mb', required=True)

	input_options = parser.add_mutually_exclusive_group(required=True)
	input_options.add_argument('-i', '--inputDirectory', type=str, metavar='/path/to/input/directory/', help='Path to directory containing the fastq files. Can be organized in separete directories by samples or all together')
	input_options.add_argument('-f', '--fastq', nargs=2, type=argparse.FileType('r'), metavar=('/path/to/input/file_1.fq.gz', '/path/to/input/file_2.fq.gz'), help='Path to Pair-End Fastq files')

	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	general_options = parser.add_argument_group('General options')
	general_options.add_argument('-o', '--outdir', type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')
	general_options.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads', required=False, default=1)
	general_options.add_argument('--json', action='store_true', help='Tells INNUca to save the results also in json format')
	general_options.add_argument('--jarMaxMemory', type=jar_max_memory, metavar='10', help='Sets the maximum RAM Gb usage by jar files (Trimmomatic and Pilon). Can also be auto or off. When auto is set, 1 Gb per thread will be used up to the free available memory', required=False, default='off')
	general_options.add_argument('--doNotUseProvidedSoftware', action='store_true', help='Tells the software to not use FastQC, Trimmomatic, SPAdes, Bowtie2, Samtools and Pilon that are provided with INNUca.py')
	# general_options.add_argument('--pairEnd_filesSeparation', nargs=2, type=str, metavar='"_left/rigth.fq.gz"', help='For unusual pair-end files separation designations, you can provide two strings containning the end of fastq files names to designate each file from a pair-end data ("_left.fq.gz" "_rigth.fq.gz" for sample_left.fq.gz sample_right.fq.gz)', required=False, default=None)
	general_options.add_argument('--keepIntermediateAssemblies', action='store_true', help='Tells INNUca to keep all the intermediate assemblies')
	general_options.add_argument('--skipEstimatedCoverage', action='store_true', help='Tells the programme to not estimate coverage depth based on number of sequenced nucleotides and expected genome size')
	general_options.add_argument('--skipTrueCoverage', action='store_true', help='Tells the programme to not run trueCoverage_ReMatCh analysis')
	general_options.add_argument('--skipFastQC', action='store_true', help='Tells the programme to not run FastQC analysis')
	general_options.add_argument('--skipTrimmomatic', action='store_true', help='Tells the programme to not run Trimmomatic')
	general_options.add_argument('--skipSPAdes', action='store_true', help='Tells the programme to not run SPAdes and consequently Pilon correction, Assembly Mapping check and MLST analysis (SPAdes contigs required)')
	general_options.add_argument('--skipAssemblyMapping', action='store_true', help='Tells the programme to not run Assembly Mapping check')
	general_options.add_argument('--skipPilon', action='store_true', help='Tells the programme to not run Pilon correction and consequently Assembly Mapping check (bam files required)')
	general_options.add_argument('--skipMLST', action='store_true', help='Tells the programme to not run MLST analysis')
	general_options.add_argument('--runPear', action='store_true', help='Tells the programme to run Pear')

	adapters_options = parser.add_mutually_exclusive_group()
	adapters_options.add_argument('--adapters', type=argparse.FileType('r'), metavar='adaptersFile.fasta', help='Fasta file containing adapters sequences to be used in FastQC and Trimmomatic', required=False)
	adapters_options.add_argument('--doNotSearchAdapters', action='store_true', help='Tells INNUca.py to not search for adapters and clip them during Trimmomatic step')

	estimated_options = parser.add_argument_group('Estimated Coverage options')
	estimated_options.add_argument('--estimatedMinimumCoverage', type=int, metavar='N', help='Minimum estimated coverage to continue INNUca pipeline', required=False, default=15)

	trueCoverage_options = parser.add_argument_group('trueCoverage_ReMatCh options')
	trueCoverage_options.add_argument('--trueConfigFile', type=argparse.FileType('r'), metavar='species.config', help='File with trueCoverage_ReMatCh settings. Some species specific config files can be found in INNUca/modules/trueCoverage_rematch/ folder. Use those files as example files. For species with config files in INNUca/modules/trueCoverage_rematch/ folder (not pre releases versions, marked with "pre."), trueCoverage_ReMatCh will run by default, unless --skipTrueCoverage is specified. Do not use together with --skipTrueCoverage option', required=False)

	fastQC_options = parser.add_argument_group('FastQC options')
	fastQC_options.add_argument('--fastQCkeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of FastQC')

	trimmomatic_options = parser.add_argument_group('Trimmomatic options')
	trimmomatic_options.add_argument('--doNotTrimCrops', action='store_true', help='Tells INNUca.py to not cut the beginning and end of reads during Trimmomatic step (unless specified with --trimCrop or --trimHeadCrop, INNUca.py will search for nucleotide content bias at both ends and will cut by there)')
	trimmomatic_options.add_argument('--trimCrop', nargs=1, type=int, metavar='N', help='Cut the specified number of bases to the end of the maximum reads length', required=False)
	trimmomatic_options.add_argument('--trimHeadCrop', nargs=1, type=int, metavar='N', help='Trimmomatic: cut the specified number of bases from the start of the reads', required=False)
	trimmomatic_options.add_argument('--trimSlidingWindow', type=str, metavar='window:meanQuality', help='Trimmomatic: perform a sliding window trimming, cutting once the average quality within the window falls below a threshold', required=False, default='5:20')
	trimmomatic_options.add_argument('--trimLeading', type=int, metavar='N', help='Trimmomatic: cut bases off the start of a read, if below a threshold quality', required=False, default=3)
	trimmomatic_options.add_argument('--trimTrailing', type=int, metavar='N', help='Trimmomatic: cut bases off the end of a read, if below a threshold quality', required=False, default=3)
	trimmomatic_options.add_argument('--trimMinLength', type=int, metavar='N', help='Trimmomatic: drop the read if it is below a specified length', required=False, default=55)
	trimmomatic_options.add_argument('--trimKeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of Trimmomatic')

	pear_options = parser.add_argument_group('Pear options')
	pear_options.add_argument('--pearKeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of Pear')
	pear_options.add_argument('--pearMinOverlap', type=int, metavar='N', help='Minimum nucleotide overlap between read pairs for Pear assembly them into only one read (default: 2/3 of maximum reads length determine using FastQC, or Trimmomatic minimum reads length if it runs, or 33 nts)', required=False)

	spades_options = parser.add_argument_group('SPAdes options')
	spades_options.add_argument('--spadesUse_3_9', action='store_true', help='Tells INNUca.py to use SPAdes v3.9.0 instead of v.3.10.1')
	spades_options.add_argument('--spadesNotUseCareful', action='store_true', help='Tells SPAdes to only perform the assembly without the --careful option')
	spades_options.add_argument('--spadesMinContigsLength', type=int, metavar='N', help='Filter SPAdes contigs for length greater or equal than this value (default: maximum reads size or 200 bp)', required=False)
	spades_options.add_argument('--spadesMaxMemory', type=int, metavar='N', help='The maximum amount of RAM Gb for SPAdes to use (default: 2 Gb per thread will be used up to the free available memory)', required=False)
	spades_options.add_argument('--spadesMinCoverageAssembly', type=spades_cov_cutoff, metavar='10', help='The minimum number of reads to consider an edge in the de Bruijn graph during the assembly. Can also be auto or off', required=False, default=2)
	spades_options.add_argument('--spadesMinKmerCovContigs', type=int, metavar='N', help='Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value', required=False, default=2)

	spades_kmers_options = parser.add_mutually_exclusive_group()
	spades_kmers_options.add_argument('--spadesKmers', nargs='+', type=int, metavar='55 77', help='Manually sets SPAdes k-mers lengths (all values must be odd, lower than 128) (default values: reads length >= 175 [55, 77, 99, 113, 127]; reads length < 175 [21, 33, 55, 67, 77])', required=False)
	spades_kmers_options.add_argument('--spadesDefaultKmers', action='store_true', help='Tells INNUca to use SPAdes default k-mers')

	assembly_mapping_options = parser.add_argument_group('Assembly Mapping options')
	assembly_mapping_options.add_argument('--assemblyMinCoverageContigs', type=int, metavar='N', help='Minimum contigs average coverage. After mapping reads back to the contigs, only keep contigs with at least this average coverage (default: 1/3 of the assembly mean coverage or 10x)', required=False)

	assembly_options = parser.add_argument_group('Assembly options')
	assembly_options.add_argument('--saveExcludedContigs', action='store_true', help='Tells INNUca.py to save excluded contigs')

	pilon_options = parser.add_argument_group('Pilon options')
	pilon_options.add_argument('--pilonKeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of Pilon')

	args = parser.parse_args()

	if args.doNotTrimCrops and (args.trimCrop or args.trimHeadCrop):
		parser.error('Cannot use --doNotTrimCrops option with --trimCrop or --trimHeadCrop')

	if args.skipTrueCoverage and args.trueConfigFile:
		parser.error('Cannot use --skipTrueCoverage option with --trueConfigFile')

	if args.spadesKmers is not None:
		for number in args.spadesKmers:
			if number % 2 == 0 or number >= 128:
				parser.error('All k-mers values must be odd integers, lower than 128')

	if len(args.speciesExpected.split(' ')) != 2:
		parser.error('Mal-formatted species name. Should be something like "Streptococcus agalactiae"')

	return args


# For parseArguments
def spades_kmers(arguments):
	kmers = sorted(map(int, arguments))
	for number in kmers:
		if number % 2 != 0 or number >= 128:
			raise argparse.ArgumentParser.error('All k-mers values must be odd integers, lower than 128')
	return kmers


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
			argparse.ArgumentParser.error('--spadesMinCoverageAssembly must be positive integer, auto or off')
	except:
		argparse.ArgumentParser.error('--spadesMinCoverageAssembly must be positive integer, auto or off')


# For parseArguments
def jar_max_memory(argument):
	string_options = ['auto', 'off']
	for option in string_options:
		if str(argument) == option:
			return str(argument)

	try:
		argument = int(argument)
		if argument > 0:
			return argument
		else:
			argparse.ArgumentParser.error('--jarMaxMemory must be positive integer, auto or off')
	except:
		argparse.ArgumentParser.error('--jarMaxMemory must be positive integer, auto or off')


def define_jar_max_memory(jar_max_memory, threads, available_memory_GB):
	GB_per_thread = 1024 / 1024.0

	maximum_memory_GB = GB_per_thread * threads

	if str(jar_max_memory) == 'off':
		return jar_max_memory
	elif str(jar_max_memory) == 'auto':
		if maximum_memory_GB > available_memory_GB:
			print 'WARNNING: the maximum memory calculated for ' + str(threads) + ' threads (' + str(round(maximum_memory_GB, 1)) + ' GB) are higher than the available memory (' + str(round(available_memory_GB, 1)) + ' GB)!'
			print 'Setting jar maximum memory to ' + str(int(round((available_memory_GB - 0.5), 0))) + ' GB'
			return int(round((available_memory_GB - 0.5), 0))
		else:
			print 'Setting jar maximum memory to ' + str(int(round(maximum_memory_GB, 0))) + ' GB'
			return int(round(maximum_memory_GB, 0))
	else:
		return jar_max_memory


def kill_subprocess_Popen(subprocess_Popen, command):
	print 'Command run out of time: ' + str(command)
	subprocess_Popen.kill()


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None, print_comand_True):
	run_successfully = False
	if not isinstance(command, basestring):
		command = ' '.join(command)
	command = shlex.split(command)

	if print_comand_True:
		print 'Running: ' + ' '.join(command)

	if shell_True:
		command = ' '.join(command)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	not_killed_by_timer = True
	if timeout_sec_None is None:
		stdout, stderr = proc.communicate()
	else:
		timer = Timer(timeout_sec_None, kill_subprocess_Popen, args=(proc, command,))
		timer.start()
		stdout, stderr = proc.communicate()
		timer.cancel()
		not_killed_by_timer = timer.isAlive()

	if proc.returncode == 0:
		run_successfully = True
	else:
		if not print_comand_True and not_killed_by_timer:
			print 'Running: ' + str(command)
		if len(stdout) > 0:
			print 'STDOUT'
			print stdout.decode("utf-8")
		if len(stderr) > 0:
			print 'STDERR'
			print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


def runTime(start_time):
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print 'Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's'
	return round(time_taken, 2)


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
		run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program, False, None, False)
		if not run_successfully:
			listMissings.append(program + ' not found in PATH.')
		else:
			print stdout.splitlines()[0]
			if programs[program][0] is None:
				print program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0]
			else:
				if program.endswith('.jar'):
					check_version = ['java', '-jar', stdout.splitlines()[0], programs[program][0]]
					programs[program].append(stdout.splitlines()[0])
				else:
					check_version = [stdout.splitlines()[0], programs[program][0]]
				run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version, False, None, False)
				if stdout == '':
					stdout = stderr
				if program == 'bunzip2':
					version_line = stdout.splitlines()[0].rsplit(',', 1)[0].split(' ')[-1]
				elif program == 'pilon-1.18.jar':
					version_line = stdout.splitlines()[0].split(' ', 3)[2]
				else:
					version_line = stdout.splitlines()[0].split(' ')[-1]
				replace_characters = ['"', 'v', 'V', '+']
				for i in replace_characters:
					version_line = version_line.replace(i, '')
				print program + ' (' + version_line + ') found'
				if programs[program][1] == '>=':
					program_found_version = version_line.split('.')
					program_version_required = programs[program][2].split('.')
					if len(program_version_required) == 3:
						if len(program_found_version) == 2:
							program_found_version.append(0)
						else:
							program_found_version[2] = program_found_version[2].split('_')[0]
					for i in range(0, len(program_version_required)):
						if int(program_found_version[i]) < int(program_version_required[i]):
							listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
				else:
					if version_line != programs[program][2]:
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
	return listMissings, programs


def setPATHvariable(args, script_path):
	path_variable = os.environ['PATH']
	script_folder = os.path.dirname(script_path)
	# Set path to use provided softwares
	if not args.doNotUseProvidedSoftware:
		fastQC = os.path.join(script_folder, 'src', 'fastqc_v0.11.5')
		trimmomatic = os.path.join(script_folder, 'src', 'Trimmomatic-0.36')
		pear = os.path.join(script_folder, 'src', 'PEAR_v0.9.10', 'bin')
		spades = os.path.join(script_folder, 'src', 'SPAdes-3.10.1-Linux', 'bin')
		if args.spadesUse_3_9:
			spades = os.path.join(script_folder, 'src', 'SPAdes-3.9.0-Linux', 'bin')
		bowtie2 = os.path.join(script_folder, 'src', 'bowtie2-2.2.9')
		samtools = os.path.join(script_folder, 'src', 'samtools-1.3.1', 'bin')
		pilon = os.path.join(script_folder, 'src', 'pilon_v1.18')

		os.environ['PATH'] = str(':'.join([fastQC, trimmomatic, pear, spades, bowtie2, samtools, pilon, path_variable]))
	print '\n' + 'PATH variable:'
	print os.environ['PATH']


# Search for fastq files in the directory
# pairEnd_filesSeparation_list can be None
def searchFastqFiles(directory, pairEnd_filesSeparation_list, listAllFilesOnly):
	filesExtensions = ['fastq.gz', 'fq.gz']
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
		sys.exit('No fastq files were found! Make sure fastq files ends with .fastq.gz or .fq.gz, and the pair-end information is either in _R1_001. or _1. format.')

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

	# Check PE files
	samples_to_remove = []
	for sample in samples:
		if len(samples[sample]) == 1:
			print 'Only one fastq file was found: ' + str(samples[sample])
			print 'Pair-End sequencing is required. This sample will be ignored'
			samples_to_remove.append(sample)

	# In Python2.7 this works
	# samples = {k: v for k, v in samples.items() if k not in samples_to_remove}
	# For Python2.6 compatibility
	temp_samples = {}
	for k, v in samples.items():
		if k not in samples_to_remove:
			temp_samples[k] = v
	samples = temp_samples
	del temp_samples

	# Create the file structure required
	for sample in samples:
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
	indir_same_outdir = False

	if inputDirectory == outdir:
		print 'Input directory and output directory are the same'
		indir_same_outdir = True

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
		sys.exit('There is no fastq files for the samples folders provided! Make sure fastq files ends with .fastq.gz or .fq.gz, and the pair-end information is either in _R1_001. or _1. format.')
	return samples, removeCreatedSamplesDirectories, indir_same_outdir


# Remove directory
def removeDirectory(directory):
	if os.path.isdir(directory):
		shutil.rmtree(directory)


# Get script version
def scriptVersionGit(version, directory, script_path):
	print 'Version ' + version

	try:
		os.chdir(os.path.dirname(script_path))
		command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, 15, False)
		print stdout
		command = ['git', 'remote', 'show', 'origin']
		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, 15, False)
		print stdout
		os.chdir(directory)
	except:
		print 'HARMLESS WARNING: git command possibly not found. The GitHub repository information will not be obtained.'


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable


def compressionType(file_to_test):
	magic_dict = {'\x1f\x8b\x08': ['gzip', 'gunzip'], '\x42\x5a\x68': ['bzip2', 'bunzip2']}

	max_len = max(len(x) for x in magic_dict)

	with open(file_to_test, 'r') as reader:
		file_start = reader.read(max_len)

	for magic, filetype in magic_dict.items():
		if file_start.startswith(magic):
			return filetype
	return None


steps = ['FastQ_Integrity', 'first_Coverage', 'trueCoverage_ReMatCh', 'first_FastQC', 'Trimmomatic', 'second_Coverage', 'second_FastQC', 'Pear', 'SPAdes', 'Pilon', 'Assembly_Mapping', 'MLST']


def sampleReportLine(run_report):
	line = []
	for step in steps:
		run_successfully = str(run_report[step][0])
		pass_qc = 'FAIL'
		if run_report[step][1] is True:
			pass_qc = 'PASS'
		elif run_report[step][1] is None:
			pass_qc = run_report[step][3]['sample']

		if step in ('first_FastQC', 'second_FastQC') and pass_qc == 'PASS' and len(run_report[step][4]) > 0:
			pass_qc = 'WARNING'
		elif step == 'SPAdes' and pass_qc == 'FAIL' and run_report['Assembly_Mapping'][1] is True:
			pass_qc = 'WARNING'

		if step in ('FastQ_Integrity', 'Pilon'):
			l = [run_successfully, run_report[step][2]]
		elif step == 'Trimmomatic':
			l = [run_successfully, pass_qc, run_report[step][2], run_report[step][4]]
		else:
			l = [run_successfully, pass_qc, run_report[step][2]]
		line.extend(l)
	return line


def start_sample_report_file(samples_report_path):
	header = ['#samples', 'samples_runSuccessfully', 'samples_passQC', 'samples_runningTime', 'samples_fileSize']
	for step in steps:
		if step == 'FastQ_Integrity':
			l = [step + '_filesOK', step + '_runningTime']
		elif step == 'Trimmomatic':
			l = [step + '_runSuccessfully', step + '_passQC', step + '_runningTime', step + '_fileSize']
		elif step == 'Pilon':
			l = [step + '_runSuccessfully', step + '_runningTime']
		else:
			l = [step + '_runSuccessfully', step + '_passQC', step + '_runningTime']
		header.extend(l)
	with open(samples_report_path, 'wt') as report:
		out = csv.writer(report, delimiter='\t')
		out.writerow(header)


def write_sample_report(samples_report_path, sample, run_successfully, pass_qc, runningTime, fileSize, run_report):
	line = [sample, run_successfully, '', runningTime, fileSize]

	line[2] = 'PASS' if pass_qc else 'FAIL'
	warning = 0
	if line[2] == 'PASS':
		if run_report['SPAdes'][1] is False:
			line[2] = 'WARNING'
			warning = 1

	line.extend(sampleReportLine(run_report))
	with open(samples_report_path, 'at') as report:
		out = csv.writer(report, delimiter='\t')
		out.writerow(line)

	return warning


def timer(function, name):
	@functools.wraps(function)
	def wrapper(*args, **kwargs):
		print('\n' + 'RUNNING {0}\n'.format(name))
		start_time = time.time()

		results = list(function(*args, **kwargs))  # guarantees return is a list to allow .insert()

		time_taken = runTime(start_time)
		print('END {0}'.format(name))

		results.insert(2, time_taken)
		return results
	return wrapper


def write_fail_report(fail_report_path, run_report):
	with open(fail_report_path, 'wt') as writer_failReport:
		failures = []
		for step in steps:
			fail_reasons = list(run_report[step][3].values())
			if fail_reasons.count(False) < len(fail_reasons):
				failures.append('#' + step)
				for key, fail_reasons in run_report[step][3].items():
					if isinstance(fail_reasons, bool) and not fail_reasons:
						continue
					else:
						failures.append('>' + str(key))
						if isinstance(fail_reasons, (list, tuple)):
							for reasons in fail_reasons:
								failures.append(str(reasons))
						else:
							failures.append(str(fail_reasons))
		writer_failReport.write('\n'.join(failures))


def parse_free_output(free_output):
	free_output = free_output.splitlines()

	dict_free_output = {}

	memory_fields = []

	counter = 0
	for line in free_output:
		if len(line) > 0:
			fields = line.split(':', 1)

			values = []
			if len(fields) == 1:
				fields = fields[0].split(' ')

				for field in fields:
					if field != '':
						values.append(field)
			elif len(fields) == 2:
				fields = fields[1].split(' ')

				for field in fields:
					if field != '':
						values.append(field)

			if counter == 0:
				memory_fields = values
			elif counter == 1:
				dict_free_output['memory'] = {}
				for i in range(0, len(memory_fields)):
					dict_free_output['memory'][memory_fields[i]] = values[i]
			elif counter == 2 and len(values) == 2:
				dict_free_output['buffers'] = {memory_fields[1]: values[0], memory_fields[2]: values[1]}
			elif counter == 3 or (counter == 2 and len(values) == 3):
				dict_free_output['swap'] = {}
				for i in range(0, 3):
					dict_free_output['swap'][memory_fields[i]] = values[i]

			counter += 1

	return dict_free_output


def get_free_memory_free(dict_free_output):
	cached = 0

	mem_cached_found = False
	for mem_type in dict_free_output['memory']:
		if 'cache' in mem_type:
			cached = int(dict_free_output['memory'][mem_type])
			mem_cached_found = True

	if not mem_cached_found:
		print 'WARNING: it was not possible to determine the cached memory!'

	return cached


def get_free_memory_os():
	free_memory_Kb = 0
	try:
		free_memory_Kb = (os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')) / 1024.0
		print 'The strict free memory available determined is ' + str(round((free_memory_Kb / (1024.0 ** 2)), 1)) + ' Gb'
	except Exception as e:
		print e
		print 'WARNING: it was not possible to determine the free memory available!'

	return free_memory_Kb


def get_free_memory():
	free_memory_Kb = 0

	command = ['free', '-k']
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, False)

	if run_successfully:
		dict_free_output = parse_free_output(stdout)

		cached = get_free_memory_free(dict_free_output)

		try:
			free_memory_Kb = int(dict_free_output['memory']['free']) + cached
			print 'For the memory use (in Kb) below, the free memory available (free + cached) is ' + str(round((free_memory_Kb / (1024.0 ** 2)), 1)) + ' Gb'
			print dict_free_output['memory']
		except:
			print 'WARNING: it was impossible to determine the free memory using the free command!'
			free_memory_Kb = get_free_memory_os()
	else:
		print 'WARNING: it was impossible to determine the free memory using the free command!'
		free_memory_Kb = get_free_memory_os()

	return free_memory_Kb


def get_cpu_information(outdir, time_str):
	with open(os.path.join(outdir, 'cpu_information.' + time_str + '.cpu.txt'), 'wt') as writer:
		command = ['cat', '/proc/cpuinfo']
		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, False)
		if run_successfully:
			writer.write(stdout)

	with open(os.path.join(outdir, 'cpu_information.' + time_str + '.slurm.txt'), 'wt') as writer:
		for environment in sorted(os.environ):
			if environment.startswith('SLURM_'):
				writer.write('#' + environment + '\n' + os.environ[environment] + '\n')


def trace_unhandled_exceptions(func):
	@functools.wraps(func)
	def wrapped_func(*args, **kwargs):
		try:
			func(*args, **kwargs)
		except:
			print 'Exception in ' + func.__name__
			traceback.print_exc()
	return wrapped_func


def get_sequence_information(fasta_file, length_extra_seq):
	sequence_dict = {}
	headers = []

	with open(fasta_file, 'rtU') as reader:
		blank_line_found = False
		sequence_counter = 0
		temp_sequence_dict = {}
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not blank_line_found:
					if line.startswith('>'):
						if len(temp_sequence_dict) > 0:
							if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
								sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
								headers.append(temp_sequence_dict.values()[0]['header'])
							else:
								print temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0'
							temp_sequence_dict = {}

						if line[1:] in headers:
							sys.exit('Found duplicated sequence headers')

						sequence_counter += 1
						temp_sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': '', 'length': 0}
					else:
						temp_sequence_dict[sequence_counter]['sequence'] += line
						temp_sequence_dict[sequence_counter]['length'] += len(line)
				else:
					sys.exit('It was found a blank line between the fasta file above line ' + line)
			else:
				blank_line_found = True

		if len(temp_sequence_dict) > 0:
			if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
				sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
				headers.append(temp_sequence_dict.values()[0]['header'])
			else:
				print temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0'

	return sequence_dict, headers
