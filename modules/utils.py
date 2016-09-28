import subprocess
import time
import sys
import argparse
import os
import shutil
import shlex
from threading import Timer
import pickle
from functools import wraps
import csv


def parseArguments(version):
	parser = argparse.ArgumentParser(prog='INNUca.py', description='INNUca - Reads Control and Assembly', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required_options = parser.add_argument_group('Required options')
	required_options.add_argument('-i', '--inputDirectory', type=str, metavar='/path/to/input/directory/', help='Path to directory containing the fastq files. Can be organized in separete directories by samples or all together', required=True)
	required_options.add_argument('-s', '--speciesExpected', type=str, metavar='"Streptococcus agalactiae"', help='Expected species name', required=True)
	required_options.add_argument('-g', '--genomeSizeExpectedMb', type=float, metavar='2.1', help='Expected genome size in Mb', required=True)

	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	general_options = parser.add_argument_group('General options')
	general_options.add_argument('-o', '--outdir', type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')
	general_options.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads', required=False, default=1)
	general_options.add_argument('--doNotUseProvidedSoftware', action='store_true', help='Tells the software to not use FastQC, Trimmomatic, SPAdes and Samtools that are provided with INNUca.py')
	general_options.add_argument('--pairEnd_filesSeparation', nargs=2, type=str, metavar='"_left/rigth.fq.gz"', help='For unusual pair-end files separation designations, you can provide two strings containning the end of fastq files names to designate each file from a pair-end data ("_left.fq.gz" "_rigth.fq.gz" for sample_left.fq.gz sample_right.fq.gz)', required=False, default=None)
	general_options.add_argument('--skipEstimatedCoverage', action='store_true', help='Tells the programme to not estimate coverage depth based on number of sequenced nucleotides and expected genome size')
	general_options.add_argument('--skipFastQC', action='store_true', help='Tells the programme to not run FastQC analysis')
	general_options.add_argument('--skipTrimmomatic', action='store_true', help='Tells the programme to not run Trimmomatic')
	general_options.add_argument('--skipSPAdes', action='store_true', help='Tells the programme to not run SPAdes and consequently MLST analysis (requires SPAdes contigs)')
	general_options.add_argument('--skipPilon', action='store_true', help='Tells the programme to not run Pilon correction')
	general_options.add_argument('--skipMLST', action='store_true', help='Tells the programme to not run MLST analysis')

	adapters_options = parser.add_mutually_exclusive_group()
	adapters_options.add_argument('--adapters', type=argparse.FileType('r'), metavar='adaptersFile.fasta', help='Fasta file containing adapters sequences to be used in FastQC and Trimmomatic', required=False)
	adapters_options.add_argument('--doNotSearchAdapters', action='store_true', help='Tells INNUca.py to not search for adapters and clip them during Trimmomatic step')

	trimmomatic_options = parser.add_argument_group('Trimmomatic options')
	trimmomatic_options.add_argument('--doNotTrimCrops', action='store_true', help='Tells INNUca.py to not cut the beginning and end of reads during Trimmomatic step (unless specified with --trimCrop or --trimHeadCrop, INNUca.py will search for nucleotide content bias at both ends and will cut by there)')
	trimmomatic_options.add_argument('--trimCrop', nargs=1, type=int, metavar='N', help='Cut the specified number of bases to the end of the maximum reads length', required=False)
	trimmomatic_options.add_argument('--trimHeadCrop', nargs=1, type=int, metavar='N', help='Trimmomatic: cut the specified number of bases from the start of the reads', required=False)
	trimmomatic_options.add_argument('--trimSlidingWindow', type=str, metavar='window:meanQuality', help='Trimmomatic: perform a sliding window trimming, cutting once the average quality within the window falls below a threshold', required=False, default='5:20')
	trimmomatic_options.add_argument('--trimLeading', type=int, metavar='N', help='Trimmomatic: cut bases off the start of a read, if below a threshold quality', required=False, default=3)
	trimmomatic_options.add_argument('--trimTrailing', type=int, metavar='N', help='Trimmomatic: cut bases off the end of a read, if below a threshold quality', required=False, default=3)
	trimmomatic_options.add_argument('--trimMinLength', type=int, metavar='N', help='Trimmomatic: drop the read if it is below a specified length', required=False, default=55)
	trimmomatic_options.add_argument('--trimKeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of Trimmomatic')

	spades_options = parser.add_argument_group('SPAdes options')
	spades_options.add_argument('--spadesNotUseCareful', action='store_true', help='Tells SPAdes to only perform the assembly without the --careful option')
	spades_options.add_argument('--spadesMinContigsLength', type=int, metavar='N', help='Filter SPAdes contigs for length greater or equal than this value', required=False, default=200)
	spades_options.add_argument('--spadesMaxMemory', type=int, metavar='N', help='The maximum amount of RAM Gb for SPAdes to use', required=False)
	spades_options.add_argument('--spadesMinCoverageAssembly', type=spades_cov_cutoff, metavar='10', help='The minimum number of reads to consider an edge in the de Bruijn graph (or path I am not sure). Can also be auto or off', required=False, default=2)
	spades_options.add_argument('--spadesMinCoverageContigs', type=int, metavar='N', help='Minimum contigs coverage. After assembly only keep contigs with reported coverage equal or above this value', required=False, default=5)

	spades_kmers_options = parser.add_mutually_exclusive_group()
	spades_kmers_options.add_argument('--spadesKmers', nargs='+', type=int, metavar='55 77', help='Manually sets SPAdes k-mers lengths (all values must be odd, lower than 128)', required=False, default=[55, 77, 99, 113, 127])
	spades_kmers_options.add_argument('--spadesDefaultKmers', action='store_true', help='Tells INNUca to use SPAdes default k-mers')

	pilon_options = parser.add_argument_group('Pilon options')
	pilon_options.add_argument('--pilonKeepFiles', action='store_true', help='Tells INNUca.py to not remove the output of Pilon')
	pilon_options.add_argument('--pilonKeepSPAdesAssembly', action='store_true', help='Tells INNUca.py to not remove the unpolished SPAdes assembly')

	args = parser.parse_args()

	if args.doNotTrimCrops and (args.trimCrop or args.trimHeadCrop):
		parser.error('Cannot use --doNotTrimCrops option with --trimCrop or --trimHeadCrop')

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


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None):
	run_successfully = False
	if isinstance(command, basestring):
		command = shlex.split(command)
	else:
		command = shlex.split(' '.join(command))

	print 'Running: ' + ' '.join(command)
	if shell_True:
		command = ' '.join(command)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if timeout_sec_None is None:
		stdout, stderr = proc.communicate()
	else:
		timer = Timer(timeout_sec_None, proc.kill)
		timer.start()
		stdout, stderr = proc.communicate()
		timer.cancel()

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
		run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program, False, None)
		if not run_successfully:
			listMissings.append(program + ' not found in PATH.')
		else:
			if programs[program][0] is None:
				print program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0]
			else:
				if program.endswith('.jar'):
					check_version = ['java', '-jar', stdout.splitlines()[0], programs[program][0]]
					programs[program].append(stdout.splitlines()[0])
				else:
					check_version = [stdout.splitlines()[0], programs[program][0]]
				run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version, False, None)
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
	return listMissings, programs


def setPATHvariable(doNotUseProvidedSoftware, script_path):
	path_variable = os.environ['PATH']
	script_folder = os.path.dirname(script_path)
	# Set path to use provided softwares
	if not doNotUseProvidedSoftware:
		fastQC = os.path.join(script_folder, 'src', 'fastqc_v0.11.5')
		trimmomatic = os.path.join(script_folder, 'src', 'Trimmomatic-0.36')
		spades = os.path.join(script_folder, 'src', 'SPAdes-3.9.0-Linux', 'bin')
		bowtie2 = os.path.join(script_folder, 'src', 'bowtie2-2.2.9')
		samtools = os.path.join(script_folder, 'src', 'samtools-1.3.1', 'bin')
		pilon = os.path.join(script_folder, 'src', 'pilon_v1.18')

		os.environ['PATH'] = str(':'.join([fastQC, trimmomatic, spades, bowtie2, samtools, pilon, path_variable]))
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
		sys.exit('There is no fastq files for the samples folders provided!')
	return samples, removeCreatedSamplesDirectories, indir_same_outdir


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


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable


# Extract compression type
def compressionType(file):
	magic_dict = {'\x1f\x8b\x08': 'gz', '\x42\x5a\x68': 'bz2'}

	max_len = max(len(x) for x in magic_dict)

	with open(file, 'r') as reader:
		file_start = reader.read(max_len)

	for magic, filetype in magic_dict.items():
		if file_start.startswith(magic):
			return filetype
	return None


steps = ('FastQ_Integrity', 'first_Coverage', 'first_FastQC', 'Trimmomatic', 'second_Coverage', 'second_FastQC', 'SPAdes', 'Pilon', 'MLST')


def sampleReportLine(run_report):
	line = []
	for step in steps:
		if step in ('FastQ_Integrity', 'Pilon'):
			l = [run_report[step][0], run_report[step][2]]
		else:
			pass_qc = 'PASS' if run_report[step][1] else 'FAIL'
			l = [run_report[step][0], pass_qc, run_report[step][2]]
		line.extend(l)
	return line


def start_sample_report_file(samples_report_path):
	header = ['#samples', 'samples_runSuccessfully', 'samples_passQC', 'samples_runningTime', 'samples_fileSize']
	for step in steps:
		if step == 'FastQ_Integrity':
			l = [step + '_filesOK', step + '_runningTime']
		elif step == 'Pilon':
			l = [step + '_runSuccessfully', step + '_runningTime']
		else:
			l = [step + '_runSuccessfully', step + '_passQC', step + '_runningTime']
		header.extend(l)
	with open(samples_report_path, 'wt') as report:
		out = csv.writer(report, delimiter='\t')
		out.writerow(header)


def write_sample_report(samples_report_path, sample, run_successfully, pass_qc, runningTime, fileSize, run_report):
	line = [sample, run_successfully, 'PASS' if pass_qc else 'FAIL', runningTime, fileSize]
	line.extend(sampleReportLine(run_report))
	with open(samples_report_path, 'at') as report:
		out = csv.writer(report, delimiter='\t')
		out.writerow(line)


def timer(function, name):
	@wraps(function)
	def wrapper(*args, **kwargs):
		print('RUNNING {0}\n'.format(name))
		start_time = time.time()

		results = list(function(*args, **kwargs))  # guarantees return is a list to allow .insert()

		time_taken = runTime(start_time)
		print('END {0}\n'.format(name))

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


def getJarPath():
	pass
