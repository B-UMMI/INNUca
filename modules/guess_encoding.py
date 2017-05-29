import utils
import os
import multiprocessing

"""
guess-encoding.py

The original guess-encoding.py script can be found in
<https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py>.

It was originally a software licensed under the MIT License reproduced bellow:

MIT License

Copyright (c) 2009-2011 Brent Pedersen, Haibao Tang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


"""
	awk 'NR % 4 == 0' your.fastq | python %prog [options]

guess the encoding of a stream of qual lines.
"""

encoding = {'Sanger': [33, (33, 73)], 'Solexa': [64, (59, 104)], 'Illumina-1.3': [64, (64, 104)], 'Illumina-1.5': [64, (66, 105)], 'Illumina-1.8': [33, (33, 74)]}


def get_qual_range(qual_str):
	vals = [ord(c) for c in qual_str]
	return min(vals), max(vals)


def get_encodings_in_range(rmin, rmax):
	valid_encodings = []
	for encoding_type, [phred, (emin, emax)] in encoding.items():
		if rmin >= emin and rmax <= emax:
			if encoding_type == 'Illumina-1.8':
				if rmax == emax:
					valid_encodings.append(['Illumina-1.8', 33])
				else:
					valid_encodings.append(['Sanger', 33])
			else:
				valid_encodings.append([encoding_type, phred])
	return valid_encodings if len(valid_encodings) > 0 else None


@utils.trace_unhandled_exceptions
def guess_encoding(fastq, number_reads_access_None_all, outdir):
	gmin, gmax = 99, 0
	valid_encodings = None
	reads_length = []
	with open(fastq, 'rtU') as reader:
		for i, line in enumerate(reader):
			if number_reads_access_None_all is None or (i + 1) / 4 <= number_reads_access_None_all:
				if (i + 1) % 4 == 0:
					if len(line) > 0:
						reads_length.append(len(line.splitlines()[0]))
						lmin, lmax = get_qual_range(line.splitlines()[0])
						if lmin < gmin or lmax > gmax:
							gmin, gmax = min(lmin, gmin), max(lmax, gmax)
							valid_encodings = get_encodings_in_range(gmin, gmax)

	utils.saveVariableToPickle([fastq, valid_encodings, min(reads_length) if len(reads_length) > 0 else None, max(reads_length) if len(reads_length) > 0 else None], outdir, 'encoding' + '.' + os.path.splitext(os.path.basename(fastq))[0])


def gather_data_together(data_directory):
	data = {}

	files = [f for f in os.listdir(data_directory) if not f.startswith('.') and os.path.isfile(os.path.join(data_directory, f))]
	for file_found in files:
		if file_found.startswith('encoding.') and file_found.endswith('.pkl'):
			file_path = os.path.join(data_directory, file_found)

			fastq, valid_encodings, min_reads_length, max_reads_length = utils.extractVariableFromPickle(file_path)
			data[fastq] = {'valid_encodings': valid_encodings, 'min_reads_length': min_reads_length, 'max_reads_length': max_reads_length}

			os.remove(file_path)

	return data


def get_final_encoding(encoding_data):
	possible_encoding = {}
	for fastq in encoding_data:
		if encoding_data[fastq] is not None and encoding_data[fastq]['valid_encodings'] is not None:
			for i in range(0, len(encoding_data[fastq]['valid_encodings'])):
				if encoding_data[fastq]['valid_encodings'][i][0] not in possible_encoding:
					possible_encoding[encoding_data[fastq]['valid_encodings'][i][0]] = 0
				possible_encoding[encoding_data[fastq]['valid_encodings'][i][0]] += 1

	final_encoding = []
	if len(possible_encoding) > 0:
		if possible_encoding.values().count(max(possible_encoding.values())) > 0:
			for encoding_type in possible_encoding:
				if possible_encoding[encoding_type] == max(possible_encoding.values()):
					final_encoding = [encoding_type, encoding[encoding_type][0]]
		else:
			final_encoding = None
	else:
		final_encoding = None

	return final_encoding


def determine_min_max_reads_length(encoding_data):
	x = [encoding_data[fastq]['min_reads_length'] for fastq in encoding_data if encoding_data[fastq]['min_reads_length'] is not None]
	min_reads_length = min(x) if len(x) > 0 else None

	x = [encoding_data[fastq]['max_reads_length'] for fastq in encoding_data if encoding_data[fastq]['max_reads_length'] is not None]
	max_reads_length = max(x) if len(x) > 0 else None

	return min_reads_length, max_reads_length


def fastq_files_enconding(fastq_files_list, number_reads_access_None_all, outdir, threads):
	pool = multiprocessing.Pool(processes=threads)
	for fastq in fastq_files_list:
		pool.apply_async(guess_encoding, args=(fastq, number_reads_access_None_all, outdir,))
	pool.close()
	pool.join()

	encoding_data = gather_data_together(outdir)

	final_encoding = get_final_encoding(encoding_data)

	min_reads_length, max_reads_length = determine_min_max_reads_length(encoding_data)

	return final_encoding, min_reads_length, max_reads_length
