#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
combine_reports.py - Combine INNUca reports
<https://github.com/miguelpmachado/manipulateFasta/>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: November 16, 2016

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details at <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import time

version = '0.1'


def get_gene_information(trueCoverage_report_file):
    gene = False

    gene_info = {}

    with open(trueCoverage_report_file, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if line.startswith('#'):
                    if line.startswith('general', 1):
                        gene = False
                    elif line.startswith('gene', 1):
                        gene = True
                    else:
                        gene = False
                else:
                    if gene:
                        line = line.split('\t')
                        gene_info[line[0]] = {'percentage_gene_coverage': float(line[1]), 'gene_mean_read_coverage': float(line[2]), 'percentage_gene_low_coverage': float(line[3]), 'number_positions_multiple_alleles': int(line[4])}

    return gene_info


def get_list_genes(gene_info):
    genes = []

    for gene in gene_info:
        genes.append(gene)

    return genes


def format_gene_info(gene_specific_info, minimum_gene_coverage):
    info = None
    if gene_specific_info['percentage_gene_coverage'] >= minimum_gene_coverage:
        if gene_specific_info['number_positions_multiple_alleles'] == 0:
            info = str(gene_specific_info['gene_mean_read_coverage'])
        else:
            info = 'multiAlleles_' + str(gene_specific_info['gene_mean_read_coverage'])
    else:
        info = 'absent_' + str(gene_specific_info['gene_mean_read_coverage'])

    return info


def combine_reports(args):
    innucaOut = os.path.abspath(args.innucaOut)
    outdir = os.path.abspath(args.outdir)
    check_create_directory(outdir)

    directories = [d for d in os.listdir(innucaOut) if not d.startswith('.') and os.path.isdir(os.path.join(innucaOut, d, ''))]

    genes = []
    results = {}

    if len(directories) == 0:
        sys.exit('No samples found')

    # Get gene list
    files = [f for f in os.listdir(os.path.join(innucaOut, directories[0], '')) if not f.startswith('.') and os.path.isfile(os.path.join(innucaOut, directories[0], f))]
    for file_found in files:
        file_path = os.path.join(innucaOut, directories[0], file_found)
        if file_found == 'trueCoverage_report.txt':
            gene_info = get_gene_information(file_path)
            genes = get_list_genes(gene_info)
            break

    if len(genes) == 0:
        sys.exit('No genes were found in one of the trueCoverage_report.txt file')

    for directory in directories:
        sample = directory
        print '\n' + sample
        directory = os.path.join(os.path.join(innucaOut, directory, ''))

        files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]
        if len(files) == 0:
            print 'No files found! Continue to the next sample'
        else:
            for file_found in files:
                file_path = os.path.join(directory, file_found)
                if file_found == 'trueCoverage_report.txt':
                    gene_info = get_gene_information(file_path)
                    if len(gene_info) != len(genes):
                        sys.exit('Different number of genes found')
                    results[sample] = {}
                    for gene in gene_info:
                        info = format_gene_info(gene_info[gene], args.minimum_gene_coverage)
                        results[sample][gene] = info
                    break

    if len(results) == 0:
        sys.exit('No results were found')

    print '\n' + 'Writing results...'
    with open(os.path.join(outdir, str('combine_trueCoverage_reports.' + time.strftime("%Y%m%d-%H%M%S") + '.tab')), 'wt') as writer:
        writer.write('#sample' + '\t' + '\t'.join(genes) + '\n')
        for sample in results:
            sample_data = [results[sample][gene] for gene in genes]
            writer.write(sample + '\t' + '\t'.join(sample_data) + '\n')

    print '\n' + 'DONE: ' + str(len(results)) + ' samples analysed.'


def check_create_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


def main():

    parser = argparse.ArgumentParser(prog='python combine_trueCoverage_reports.py', description="Combine trueCoverage_ReMatCh module reports in respect to gene information.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--innucaOut', type=str, metavar='/path/to/INNUca/output/directory/', help='Path to INNUca output directory', required=True)

    parser_optional = parser.add_argument_group('Facultative options')
    parser_optional.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/', help='Path to where to store the outputs', required=False, default='.')
    parser_optional.add_argument('--minimum_gene_coverage', type=int, metavar='80', help='Minimum percentage of sequence length (with a minimum of read depth to consider a position to be present) to determine whether a gene is present.', required=False, default=80)

    parser.set_defaults(func=combine_reports)

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
