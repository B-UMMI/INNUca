#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
kraken.py - Runs Kraken on fasta and fastq files, parse the results and evaluate the outputs
<https://github.com/B-UMMI/INNUca/modules/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: September 20, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os.path
from utils import runCommandPopenCommunicate as utils_run_command
from utils import runTime as utils_run_time
from utils import checkPrograms as utils_check_programs
from utils import timer as utils_timer
from functools import partial
import argparse
import time

version = '1.0'


def find_compression_type(file_to_test):
    """
    Determines the compression type of a file

    Parameters
    ----------
    file_to_test : str
        Path fot the file from which the compression type will be determined

    Returns
    -------
    compression_type : str or None
        String with the compression type. If none of the compression types searched is found, None is returned
    """
    magic_dict = {b'\x1f\x8b\x08': 'gzip',
                  b'\x42\x5a\x68': 'bzip2',
                  b'\x50\x4b\x03\x04': 'zip',
                  b'\xfd\x37\x7a\x58\x5a\x00': 'xz',
                  b'\x37\x7a\xbc\xaf\x27\x1c': '7zip'}

    max_len = max(len(x) for x in magic_dict)

    with open(file_to_test, 'rb') as reader:
        file_start = reader.read(max_len)

    compression_type = None
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            compression_type = filetype

    return compression_type


def get_kraken_version():
    """
    Determine Kraken version. If kraken2 exists, it will be the default.

    Returns
    -------
    version_kraken : int or None
        1 if only first Kraken (or zero version) was found, 2 if kraken2 was found, or None if none of those was found
    """
    version_kraken = None

    command = ['which', 'kraken']
    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=False)
    if run_successfully:
        version_kraken = 1

    command[1] = 'kraken2'
    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=False)
    if run_successfully:
        version_kraken = 2

    return version_kraken


def kraken_compression_type(file_to_test):
    """
    Get compression type and check if compatible with Kraken (gzip, bzip2 or uncompressed allowed)

    Parameters
    ----------
    file_to_test : str
        Path fot the file from which the compression type will be determined

    Returns
    -------
    compression_type : str
        String with the compression type. If none of the compression types searched is found, None is returned
    """
    compression_allowed = ['gzip', 'bzip2', None]
    compression_type = find_compression_type(file_to_test)
    if compression_type not in compression_allowed:
        sys.exit('Compression type found ({compression_type}) not allowed.\n'
                 'Use only: {compression_allowed}'.format(compression_type=compression_type,
                                                          compression_allowed=compression_allowed))
    return compression_type


def run_kraken_main(files_to_classify, kraken_db, files_type, outdir, version_kraken, db_mem=False, quick=False,
                    min_base_quality=10, threads=1):
    """
    Run Kraken for data classification

    Parameters
    ----------
    files_to_classify : list
        List with files to be classified by Kraken. Can be one fasta or up to two fastq files
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    files_type : str
        Type of the files to be classified: fasta or fastq
    outdir : str
        Path to the output directory
    version_kraken : int or None
        1 if only first Kraken (or zero version) was found, 2 if kraken2 was found, or None if none of those was found
    db_mem : bool, default False
        True if want to load the Kraken DB into memory before run, else False
    quick : bool, default False
        True if want to do a quick operation and only use the first hits
    min_base_quality : int, default 10
        Minimum base quality used in classification. Only used with fastq files and kraken2.
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if Kraken ran successfully or not
    kraken_output : str or None
        Path to Kraken output file. If running kraken2, None is returned
    kraken_report : str or None
        Path to Kraken report (results) file in case of kraken2, else None
    """
    files_type_options = ['fastq', 'fasta']
    if files_type not in files_type_options:
        raise ValueError("Invalid files type. Expected one of: %s" % files_type_options)

    kraken_output = os.path.join(outdir, 'kraken.{db}.out'.format(db=os.path.basename(kraken_db)))
    kraken_report = None

    command = ['kraken', '', '', '--db', kraken_db, '--threads', str(threads), '--output', kraken_output, '',
               '--{type}-input'.format(type=files_type), '', '', '', '', '', '', '', ' '.join(files_to_classify)]

    if version_kraken == 2:
        command[0] = 'kraken2'
        command[8] = '-'
        command[10] = ''
        kraken_output = None
        kraken_report = os.path.join(outdir, 'kraken_report.{db}.txt'.format(db=os.path.basename(kraken_db)))
        command[11] = '--report'
        command[12] = kraken_report
        """
        Didn't get what this confidence mean
          --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                                  in [0, 1].

        At present, we have not yet developed a confidence score with a probabilistic interpretation for Kraken 2.
         However, we have developed a simple scoring scheme that has yielded good results for us, and we've made that
          available in Kraken 2 through use of the --confidence option to kraken2. The approach we use allows a user to
           specify a threshold score in the [0,1] interval; the classifier then will adjust labels up the tree until the
            label's score (described below) meets or exceeds that threshold. If a label at the root of the taxonomic
             tree would not have a score exceeding the threshold, the sequence is called unclassified by Kraken 2 when
              this threshold is applied.

        A sequence label's score is a fraction C/Q, where C is the number of k-mers mapped to LCA values in the clade
         rooted at the label, and Q is the number of k-mers in the sequence that lack an ambiguous nucleotide (i.e.,
          they were queried against the database). Consider the example of the LCA mappings in Kraken 2's output given
           earlier:

        "562:13 561:4 A:31 0:1 562:3" would indicate that:
        
        the first 13 k-mers mapped to taxonomy ID #562
        the next 4 k-mers mapped to taxonomy ID #561
        the next 31 k-mers contained an ambiguous nucleotide
        the next k-mer was not in the database
        the last 3 k-mers mapped to taxonomy ID #562

        In this case, ID #561 is the parent node of #562. Here, a label of #562 for this sequence would have a score of
         C/Q = (13+3)/(13+4+1+3) = 16/21. A label of #561 would have a score of C/Q = (13+4+3)/(13+4+1+3) = 20/21. If a
          user specified a --confidence threshold over 16/21, the classifier would adjust the original label from #562
           to #561; if the threshold was greater than 20/21, the sequence would become unclassified.
        """
        # command[13] = '--confidence'
        # command[14] = '1'
        if files_type == 'fastq':
            command[15] = '--minimum-base-quality'
            command[16] = str(min_base_quality)

    if len(files_to_classify) == 0:
        sys.exit('No files provided for classification.')
    elif len(files_to_classify) <= 2:
        if files_type == 'fastq' and len(files_to_classify) == 2:
            command[17] = '--paired'
        elif files_type == 'fasta':
            if len(files_to_classify) == 2:
                sys.exit('{n} files provided for classification. Maximum of 1 file for fasta is'
                         ' allowed.'.format(n=len(files_to_classify)))
    elif len(files_to_classify) > 2:
        sys.exit('{n} files provided for classification. Maximum of 2 files for fastq or 1 file for fasta are'
                 ' allowed.'.format(n=len(files_to_classify)))

    compression_type = kraken_compression_type(files_to_classify[0])
    if compression_type is not None:
        command[9] = '--{type}-compressed'.format(type=compression_type)
    del compression_type

    if quick:
        command[1] = '--quick'

    if db_mem and version_kraken == 1:
        command[2] = '--preload'
    elif not db_mem and version_kraken == 2:
        command[2] = '--memory-mapping'

    run_successfully, _, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                               print_comand_True=True)

    return run_successfully, kraken_output, kraken_report


def run_kraken_report(kraken_db, kraken_output, outdir):
    """
    Get the Kraken report from kraken run

    Parameters
    ----------
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    kraken_output : str
        Path to Kraken output file
    outdir : str
        Path to the output directory

    Returns
    -------
    run_successfully : bool
        Boolean stating if Kraken ran successfully or not
    kraken_results : str
        String with Kraken report
    """

    command = ['kraken-report', '--db', kraken_db, kraken_output]
    run_successfully, kraken_results, _ = utils_run_command(command=command, shell_True=False, timeout_sec_None=None,
                                                            print_comand_True=True)
    if run_successfully:
        with open(os.path.join(outdir, 'kraken_report.{db}.txt'.format(db=os.path.basename(kraken_db))),
                  'wt') as writer:
            writer.write(kraken_results)

    return run_successfully, kraken_results


def read_kraken_report(kraken_report):
    """
    Read Kraken report into a variable

    Parameters
    ----------
    kraken_report : str
        Path to Kraken report (results) file

    Returns
    -------
    kraken_results : str
        String with Kraken report
    """
    kraken_results = None

    with open(kraken_report, 'rt') as reader:
        kraken_results = reader.read()

    return kraken_results


def parse_kraken_results(kraken_results):
    """
    Parse the Kraken report

    Parameters
    ----------
    kraken_results : str
        String with Kraken report. If using a file to store the report, save the file into a variable

    Returns
    -------
    results_parsed : dict
        Dictionary with the percentage of fragments covered by the clade rooted at this taxon. species classifications
        are stored (or alternatively the information for genus is stored if no species for a given genus exists). Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    unclassified : float
        Percentage of unclassified fragments
    """

    results_parsed = {}
    genus = {}
    species = {}
    unclassified = 0
    for line in kraken_results.splitlines():
        if len(line) > 0:
            fields = line.split('\t')
            if fields[3].strip() == 'U':
                unclassified = float(fields[0].strip())
            elif fields[3].strip() == 'G':
                if len(genus) == 0:
                    genus[('genus', fields[4].strip(), fields[5].strip())] = float(fields[0].strip())
                else:
                    if len(species) == 0:
                        results_parsed.update(genus)
                    else:
                        results_parsed.update(species)
                    genus = {('genus', fields[4].strip(), fields[5].strip()): float(fields[0].strip())}
                    species = {}
            elif fields[3].strip() == 'S':
                species[('species', fields[4].strip(), fields[5].strip())] = float(fields[0].strip())
            elif not fields[3].strip().startswith(('G', 'S')):
                if len(species) == 0 and len(genus) > 0:
                    results_parsed.update(genus)
                elif len(species) > 0:
                    results_parsed.update(species)
                genus = {}
                species = {}
    if len(species) == 0:
        results_parsed.update(genus)
    else:
        results_parsed.update(species)

    return results_parsed, unclassified


def clean_kraken_results(results_parsed, min_percent_covered=None):
    """
    Clean the Kraken results

    Parameters
    ----------
    results_parsed : dict
        Dictionary with the percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    min_percent_covered : float or None, default None
        Minimum percentage of fragments covered to consider the taxon. If None is provided, the "maximum percentage of
        fragments covered found in results_parsed" / 100 will be used.

    Returns
    -------
    cleaned_results : dict
        Dictionary with the cleaned percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    """
    # for k, v in results_parsed.items():
    #     results_parsed[k] = float(v)

    cleaned_results = {}
    if len(results_parsed) > 0:
        max_percent_covered = sorted(results_parsed.items(), key=lambda x: x[1], reverse=True)[0][1]

        if min_percent_covered is None:
            min_percent_covered = max_percent_covered / 100

        del max_percent_covered

        for species_info, percent_covered in results_parsed.items():
            if percent_covered >= min_percent_covered:
                cleaned_results[species_info] = percent_covered

    return cleaned_results


def run_kraken_module(files_to_classify, kraken_db, files_type, outdir, version_kraken, db_mem=False, quick=False,
                      min_percent_covered=None, min_base_quality=10, threads=1):
    """
    Run Kraken, parse and clean the results

    Parameters
    ----------
    files_to_classify : list
        List with files to be classified by Kraken. Can be one fasta or up to two fastq files
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    files_type : str
        Type of the files to be classified: fasta or fastq
    outdir : str
        Path to the output directory
    version_kraken : int or None
        1 if only first Kraken (or zero version) was found, 2 if kraken2 was found, or None if none of those was found
    db_mem : bool, default False
        True if want to load the Kraken DB into memory before run, else False
    quick : bool, default False
        True if want to do a quick operation and only use the first hits
    min_percent_covered : float or None, default None
        Minimum percentage of fragments covered to consider the taxon. If None is provided, the ["maximum percentage of
        fragments covered found (excluding unclassified category)" / 100] will be used.
    min_base_quality : int, default 10
        Minimum base quality used in classification. Only used with fastq files and kraken2.
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if Kraken ran successfully or not
    cleaned_results : dict
        Dictionary with the cleaned percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    unclassified : float
        Percentage of unclassified fragments
    """

    cleaned_results = {}
    unclassified = 0

    run_successfully, kraken_output, kraken_report = run_kraken_main(files_to_classify=files_to_classify,
                                                                     kraken_db=kraken_db, files_type=files_type,
                                                                     outdir=outdir, version_kraken=version_kraken,
                                                                     db_mem=db_mem, quick=quick,
                                                                     min_base_quality=min_base_quality, threads=threads)
    kraken_results = None
    if run_successfully:
        if version_kraken == 2:
            kraken_results = read_kraken_report(kraken_report)
        else:
            run_successfully, kraken_results = run_kraken_report(kraken_db=kraken_db, kraken_output=kraken_output,
                                                                 outdir=outdir)

    if kraken_output is not None and os.path.isfile(kraken_output):
        os.remove(kraken_output)
    del kraken_output

    if run_successfully:
        results_parsed, unclassified = parse_kraken_results(kraken_results=kraken_results)

        if len(results_parsed) > 0:
            cleaned_results = clean_kraken_results(results_parsed=results_parsed,
                                                   min_percent_covered=min_percent_covered)
        del results_parsed
    del kraken_results

    return run_successfully, cleaned_results, unclassified


def get_most_abundant_taxon(cleaned_results):
    """
    Get taxonomy information and percentage of fragments assigned to the most abundant taxon

    Parameters
    ----------
    cleaned_results : dict
        Dictionary with the cleaned percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage

    Returns
    -------
    taxonomy_level : str or None
        Taxonomy level: 'genus', 'species' or None
    taxonomy_id : str or None
        Taxonomy ID
    taxon_name : str or None
        Species name or genus name
    taxon_percentage : float or None
        Percentage of fragments assigned to this taxon
    """

    taxon_level = None
    taxon_id = None
    taxon_name = None
    taxon_percentage = None
    if len(cleaned_results) > 0:
        sorted_results = sorted(cleaned_results.items(), key=lambda x: x[1], reverse=True)
        taxon_level = sorted_results[0][0][0]
        taxon_id = sorted_results[0][0][1]
        taxon_name = sorted_results[0][0][2]
        taxon_percentage = sorted_results[0][1]

    return taxon_level, taxon_id, taxon_name, taxon_percentage


def evaluate_multispecies_unknown(cleaned_results, unclassified, max_unclassified_frag=None):
    """

    Parameters
    ----------
    cleaned_results : dict
        Dictionary with the cleaned percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    unclassified : float
        Percentage of unclassified fragments
    max_unclassified_frag : float or None, default None
        Maximum percentage of unclassified fragments. If None is provided, the [100 - "maximum percentage of
        fragments covered found (excluding unclassified category)" / 10] will be used.

    Returns
    -------
    multispecies : bool
        States if multiples species was found. If more than one taxon is found in cleaned_results, True is returned.
    unknown_fragments : bool
        States if there are too many unknown (unclassified) fragments
    """

    multispecies = False
    unknown_fragments = False

    if len(cleaned_results) > 1:
        multispecies = True

    if len(cleaned_results) > 0:
        if max_unclassified_frag is None:
            max_unclassified_frag = (100 - sorted(cleaned_results.items(), key=lambda x: x[1], reverse=True)[0][1]) / 10

        if unclassified > max_unclassified_frag:
            unknown_fragments = True
    else:
        unknown_fragments = True

    return multispecies, unknown_fragments


def write_reports(multispecies, unknown_fragments, cleaned_results, unclassified, kraken_db, outdir,
                  most_abundant_taxon=None):
    """
    Write Kraken results

    Parameters
    ----------
    multispecies : bool
        States if multiples species was found
    unknown_fragments : bool
        States if there are too many unknown (unclassified) fragments
    cleaned_results : dict
        Dictionary with the cleaned percentage of fragments covered by the clade rooted at this taxon. Keys
        are tuples of (taxonomy_id, species_name) and values are floats with the percentage
    unclassified : float
        Percentage of unclassified fragments
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    outdir : str
        Path to the output directory
    most_abundant_taxon : tuple or None, default None
        Tuple with three strings for taxon level, id and name (in this order), and a fourth float element with the
        percentage of fragments assigned to this taxon. If None is provided, those three
        information will be determined.

    Returns
    -------

    """
    results_classification = os.path.join(outdir,
                                          'kraken_results.classification.{}.tab'.format(os.path.basename(kraken_db)))
    # with open(results_classification, 'wt', newline='\n') as writer:
    with open(results_classification, 'wt') as writer:
        writer.write('\t'.join(['#taxon_level', 'taxon_id', 'taxon_name', 'percentage_fragments_covered']) + '\n')
        writer.write('\t'.join(['unclassified', '', '', str(unclassified)]) + '\n')
        for species_info, percent_covered in cleaned_results.items():
            writer.write('\t'.join([species_info[0], species_info[1], species_info[2], str(percent_covered)]) + '\n')

    results_evaluation = os.path.join(outdir, 'kraken_results.evaluation.{}.tab'.format(os.path.basename(kraken_db)))
    # with open(results_evaluation, 'wt', newline='\n') as writer:
    with open(results_evaluation, 'wt') as writer:
        if most_abundant_taxon is None:
            most_abundant_taxon_level, most_abundant_taxon_id, most_abundant_taxon_name, most_abundant_taxon_percent = \
                get_most_abundant_taxon(cleaned_results)
        else:
            most_abundant_taxon_level = most_abundant_taxon[0]
            most_abundant_taxon_id = most_abundant_taxon[1]
            most_abundant_taxon_name = most_abundant_taxon[2]
            most_abundant_taxon_percent = most_abundant_taxon[3]

        writer.write('\t'.join(['#multiple_taxon_found', 'number_taxon_found',
                                'excessive_unknown_fragments', 'percentage_unknown_fragments',
                                'most_abundant_taxon_level', 'most_abundant_taxon_id', 'most_abundant_taxon_name',
                                'percentage_most_abundant_taxon_fragments']) + '\n')
        writer.write('\t'.join([str(multispecies), str(len(cleaned_results)),
                                str(unknown_fragments), str(unclassified),
                                most_abundant_taxon_level, most_abundant_taxon_id, most_abundant_taxon_name,
                                str(most_abundant_taxon_percent)]) + '\n')


def print_results(multispecies, unknown_fragments, most_abundant_taxon):
    """
    Print general Kraken results

    Parameters
    ----------
    multispecies : bool
        States if multiples species was found
    unknown_fragments : bool
        States if there are too many unknown (unclassified) fragments
    most_abundant_taxon : tuple
        Tuple with three strings for taxon level, id and name (in this order), and a fourth float element with the
        percentage of fragments assigned to this taxon

    Returns
    -------

    """
    print('\n'
          '--- RESULTS ---\n')

    if multispecies:
        print('Multiple species was found!')
    if unknown_fragments:
        print('Too many unknown (unclassified) fragments were found!')
    print('\n'
          'The most abundant taxon found was: {species} ({percent}% of'
          ' fragments_assigned) ({level} id {id}).\n'.format(species=most_abundant_taxon[2],
                                                             level=most_abundant_taxon[0],
                                                             id=most_abundant_taxon[1],
                                                             percent=most_abundant_taxon[3]))


def run_module(files_to_classify, kraken_db, files_type, outdir, version_kraken, db_mem=False, quick=False,
               min_percent_covered=None, max_unclassified_frag=None, min_base_quality=10, threads=1):
    """
    Runs Kraken, parse the results and evaluate the outputs

    Parameters
    ----------
    files_to_classify : list
        List with files to be classified by Kraken. Can be one fasta or up to two fastq files
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    files_type : str
        Type of the files to be classified: fasta or fastq
    outdir : str
        Path to the output directory
    version_kraken : int or None
        1 if only first Kraken (or zero version) was found, 2 if kraken2 was found, or None if none of those was found
    db_mem : bool, default False
        True if want to load the Kraken DB into memory before run, else False
    quick : bool, default False
        True if want to do a quick operation and only use the first hits
    min_percent_covered : float or None, default None
        Minimum percentage of fragments covered to consider the taxon. If None is provided, the ["maximum percentage of
        fragments covered found (excluding unclassified category)" / 100] will be used.
    max_unclassified_frag : float or None, default None
        Maximum percentage of unclassified fragments. If None is provided, the [100 - "maximum percentage of
        fragments covered found (excluding unclassified category)" / 10] will be used.
    min_base_quality : int, default 10
        Minimum base quality used in classification. Only used with fastq files and kraken2.
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if Kraken ran successfully or not
    multispecies : bool
        States if multiples species was found
    unknown_fragments : bool
        States if there are too many unknown (unclassified) fragments
    most_abundant_taxon : tuple
        Tuple with three strings for taxon level, id and name (in this order), and a fourth float element with the
        percentage of fragments assigned to this taxon
    """

    run_successfully, cleaned_results, unclassified = run_kraken_module(files_to_classify=files_to_classify,
                                                                        kraken_db=kraken_db, files_type=files_type,
                                                                        outdir=outdir, version_kraken=version_kraken,
                                                                        db_mem=db_mem, quick=quick,
                                                                        min_percent_covered=min_percent_covered,
                                                                        min_base_quality=min_base_quality,
                                                                        threads=threads)

    multispecies = False
    unknown_fragments = False
    most_abundant_taxon = (None, None, None, None)
    if run_successfully:
        multispecies, unknown_fragments = evaluate_multispecies_unknown(cleaned_results=cleaned_results,
                                                                        unclassified=unclassified,
                                                                        max_unclassified_frag=max_unclassified_frag)

        most_abundant_taxon_level, most_abundant_taxon_id, most_abundant_taxon_name, most_abundant_taxon_percent = \
            get_most_abundant_taxon(cleaned_results)
        most_abundant_taxon = (most_abundant_taxon_level, most_abundant_taxon_id, most_abundant_taxon_name,
                               most_abundant_taxon_percent)

        print_results(multispecies=multispecies, unknown_fragments=unknown_fragments,
                      most_abundant_taxon=most_abundant_taxon)

        write_reports(multispecies=multispecies, unknown_fragments=unknown_fragments, cleaned_results=cleaned_results,
                      unclassified=unclassified, kraken_db=kraken_db, outdir=outdir,
                      most_abundant_taxon=most_abundant_taxon)

    return run_successfully, multispecies, unknown_fragments, most_abundant_taxon


def rename_output_innuca(kraken_db, files_type, outdir):
    """
    Rename the output files to contain the information if it was ran in the reads or in the assembly

    Parameters
    ----------
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    files_type : str
        Type of the files to be classified: fasta or fastq
    outdir : str
        Path to the output directory

    Returns
    -------
    """

    # kraken_report
    old = os.path.join(outdir, 'kraken_report.{db}.txt'.format(db=os.path.basename(kraken_db)))
    new = os.path.join(outdir, 'kraken_report.{db}.{type}.txt'.format(db=os.path.basename(kraken_db), type=files_type))
    if os.path.isfile(old):
        os.rename(old, new)

    # kraken classification cleaned
    old = os.path.join(outdir, 'kraken_results.classification.{db}.tab'.format(db=os.path.basename(kraken_db)))
    new = os.path.join(outdir, 'kraken_results.classification.{db}.{type}.tab'.format(db=os.path.basename(kraken_db),
                                                                                      type=files_type))
    if os.path.isfile(old):
        os.rename(old, new)

    # kraken evaluation
    old = os.path.join(outdir, 'kraken_results.evaluation.{db}.tab'.format(db=os.path.basename(kraken_db)))
    new = os.path.join(outdir, 'kraken_results.evaluation.{db}.{type}.tab'.format(db=os.path.basename(kraken_db),
                                                                                  type=files_type))
    if os.path.isfile(old):
        os.rename(old, new)


kraken_timer = partial(utils_timer, name='Kraken')


@kraken_timer
def run_for_innuca(species, files_to_classify, kraken_db, files_type, outdir, version_kraken, db_mem=False, quick=False,
                   min_percent_covered=None, max_unclassified_frag=None, min_base_quality=10, threads=1):
    """
    Runs Kraken for INNUca and QA/QC the results

    Parameters
    ----------
    species : str
        Species name, e.g. "streptococcus agalactiae"
    files_to_classify : list
        List with files to be classified by Kraken. Fastq files required
    kraken_db : str
        Kraken DB name or path to the directory containing the Kraken DB
    files_type : str
        Type of the files to be classified: fasta or fastq
    outdir : str
        Path to the output directory
    version_kraken : int or None
        1 if only first Kraken (or zero version) was found, 2 if kraken2 was found, or None if none of those was found
    db_mem : bool, default False
        True if want to load the Kraken DB into memory before run, else False
    quick : bool, default False
        True if want to do a quick operation and only use the first hits
    min_percent_covered : float or None, default None
        Minimum percentage of fragments covered to consider the taxon. If None is provided, the ["maximum percentage of
        fragments covered found (excluding unclassified category)" / 100] will be used.
    max_unclassified_frag : float or None, default None
        Maximum percentage of unclassified fragments. If None is provided, the [100 - "maximum percentage of
        fragments covered found (excluding unclassified category)" / 10] will be used.
    min_base_quality : int, default 10
        Minimum base quality used in classification. Only used with fastq files and kraken2.
    threads : int, default 1
        Number of threads to be used

    Returns
    -------
    run_successfully : bool
        Boolean stating if INNUca Kraken module ran successfully or not
    pass_qc : bool
        Boolean stating if sample pass QA/QC or not
    time_taken : float
        Seconds that run_for_innuca took to run
    failing : dict
        Dictionary with the failing reasons. If sample did not fail, it is only {'sample': False}. If it failed, keys
        will be the level of failing, and values list of strings
    warnings : dict
        Dictionary with the warning reasons. If no warnings were raised, it is empty. If warnings were raised, keys
        will be the level of warnings, and values list of strings
    most_abundant_taxon_percent : float or None
        Percentage of fragments assigned to the most abundant taxon. If no taxon was found, None is returned
    """

    pass_qc = False
    failing = {}
    warnings = {}

    if os.path.isdir(kraken_db):
        kraken_db = os.path.abspath(kraken_db)
        if kraken_db.endswith('/'):
            kraken_db = kraken_db[:-1]

    species = species.lower().split(' ')

    run_successfully, multispecies, unknown_fragments, most_abundant_taxon = \
        run_module(files_to_classify=files_to_classify, kraken_db=kraken_db, files_type=files_type, outdir=outdir,
                   version_kraken=version_kraken, db_mem=db_mem, quick=quick, min_percent_covered=min_percent_covered,
                   max_unclassified_frag=max_unclassified_frag, min_base_quality=min_base_quality, threads=threads)

    rename_output_innuca(kraken_db=kraken_db, files_type=files_type, outdir=outdir)

    # QA/QC assessment
    if run_successfully:
        if most_abundant_taxon is not None:
            taxon_name = most_abundant_taxon[2].lower().split(' ')
            if most_abundant_taxon[0] == 'genus':
                if species[0] == taxon_name[0]:
                    warnings['taxon'] = ['Only genus results were found by Kraken']
                else:
                    failing['taxon'] = ['The most abundant taxon found by'
                                        ' Kraken ({kraken_genus} genus) does not match the provided'
                                        ' genus ({species_genus})'.format(kraken_genus=taxon_name[0],
                                                                          species_genus=species[0])]
            else:
                if species != taxon_name[:2]:
                    failing['taxon'] = ['The most abundant species found by'
                                        ' Kraken ({kraken_species}) does not match the provided'
                                        ' species ({species})'.format(kraken_species=' '.join(taxon_name[:2]),
                                                                      species=' '.join(species))]
        else:
            warnings['taxon'] = ['No taxon was found by Kraken']

        if unknown_fragments:
            warnings['unknown'] = ['Too many unknown (unclassified) fragments']

        if multispecies:
            warnings['multispecies'] = ['Multiple taxon found']
    else:
        failing['sample'] = ['Did not run']

    if len(failing) == 0:
        failing = {'sample': False}
        pass_qc = True
    else:
        print('Failing: {}'.format(failing))

    most_abundant_taxon_percent = most_abundant_taxon[3]

    return run_successfully, pass_qc, failing, warnings, most_abundant_taxon_percent


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 kraken.py"')

    parser = argparse.ArgumentParser(prog='kraken.py',
                                     description='Runs Kraken on fasta and fastq files, parse the results and evaluate'
                                                 ' the outputs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-f', '--files', nargs='+', type=argparse.FileType('r'), required=True,
                                 metavar='/path/to/files/to/classify.fastq.gz',
                                 help='Path to the files to be classified by Kraken. Can be one fasta or up to two'
                                      ' fastq files')
    parser_required.add_argument('-k', '--kraken', type=str, metavar='minikraken_20171013_4GB',
                                 help='Name of Kraken DB found in path, or complete path to the directory'
                                      ' containing the Kraken DB files (for example'
                                      ' /path/to/directory/minikraken_20171013_4GB)',
                                 required=True)
    parser_required.add_argument('-t', '--type', choices=['fasta', 'fastq'], type=str, metavar='fasta', required=True,
                                 help='Input files type (available options: %(choices)s)')

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the sequences will be stored (default: ./)',
                                         required=False, default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N',
                                         help='Number of threads to use (default: 1)', required=False, default=1)

    parser_optional_kraken = parser.add_argument_group('Kraken facultative options')
    parser_optional_kraken.add_argument('-m', '--memory', action='store_true',
                                        help='Set Kraken to load the DB into the memory before run')
    parser_optional_kraken.add_argument('-q', '--quick', action='store_true',
                                        help='Set Kraken to do a quick operation and only use the first hits')
    parser_optional_kraken.add_argument('--min_cov', type=float, metavar='1.5', required=False,
                                        help='Minimum percentage of fragments covered to consider the taxon. If nothing'
                                             ' is specified, the hundredth of the taxon found (species or genus if no'
                                             ' species are available for a given genus) with higher percentage of'
                                             ' fragments covered (excluding unclassified category) will be used.')
    parser_optional_kraken.add_argument('--max_unclass', type=float, metavar='1.5', required=False,
                                        help='Maximum percentage of unclassified fragments allowed. If nothing is'
                                             ' specified, the tenth of 100 minus the percentage of fragments of the'
                                             ' taxon found (species or genus if no  species are available for a given'
                                             ' genus) with higher percentage of fragments covered (excluding'
                                             ' unclassified category) will be used.')
    parser_optional_kraken.add_argument('--min_quality', type=int, metavar='N',
                                        help='Using fastq files, sets the minimum base quality to be used in'
                                             ' classification (default: 10)', required=False, default=10)

    args = parser.parse_args()

    msg = []
    if args.min_cov is not None:
        if args.min_cov < 0 or args.min_cov > 100:
            msg.append('--min_cov should be a value between [0, 100]')
    if args.max_unclass is not None:
        if args.max_unclass < 0 or args.max_unclass > 100:
            msg.append('--max_unclass should be a value between [0, 100]')

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

    start_time = time.time()

    print('\n'
          '===>  RUNNING  kraken.py  <===\n')

    version_kraken = get_kraken_version()
    if version_kraken == 2:
        missing_programs, _ = utils_check_programs({'kraken2': ['--version', '>=', '2.0.6']})
    elif version_kraken == 1:
        missing_programs, _ = utils_check_programs({'kraken': ['--version', '>=', '0.10.6'],
                                                    'kraken-report': ['--version', '>=', '0.10.6']})
    else:
        print()
        sys.exit('Kraken was not found in PATH')

    print()

    if len(missing_programs) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missing_programs))

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    args.files = [os.path.abspath(file_classify.name) for file_classify in args.files]

    if os.path.isdir(args.kraken):
        args.kraken = os.path.abspath(args.kraken)
        if args.kraken.endswith('/'):
            args.kraken = args.kraken[:-1]

    run_successfully, _, _, _ = run_module(files_to_classify=args.files, kraken_db=args.kraken,
                                           files_type=args.type, outdir=args.outdir, version_kraken=version_kraken,
                                           db_mem=args.memory, quick=args.quick, min_percent_covered=args.min_cov,
                                           max_unclassified_frag=args.max_unclass, min_base_quality=args.min_quality,
                                           threads=args.threads)

    _ = utils_run_time(start_time)

    if not run_successfully:
        sys.exit('Something went wrong!')


if __name__ == "__main__":
    main()
