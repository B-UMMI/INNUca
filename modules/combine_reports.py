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

version = '0.3'


def combine_reports(innucaOut, outdir, json, time_str, number_samples):
    if time_str is None:
        time_str = time.strftime("%Y%m%d-%H%M%S")

    innucaOut = os.path.abspath(innucaOut)
    outdir = os.path.abspath(outdir)
    check_create_directory(outdir)

    files = [f for f in os.listdir(innucaOut) if not f.startswith('.') and os.path.isfile(os.path.join(innucaOut, f))]
    samples_report = []
    for file_found in files:
        file_found = os.path.join(innucaOut, file_found)
        if os.path.basename(file_found).startswith('samples_report'):
            samples_report.append(file_found)
    if len(samples_report) == 0:
        samples_report = None
    elif len(samples_report) == 1:
        samples_report = samples_report[0]
    else:
        samples_report = sorted(samples_report, reverse=True)[0]

    directories = [d for d in os.listdir(innucaOut) if not d.startswith('.') and os.path.isdir(os.path.join(innucaOut, d, ''))]

    results = {}

    if len(directories) == 0:
        sys.exit('No samples found')

    fields = ['#samples', 'first_coverage', 'trueCoverage_absent_genes', 'trueCoverage_multiple_alleles', 'trueCoverage_sample_coverage', 'second_Coverage', 'pear_assembled_reads', 'pear_unassembled_reads', 'pear_dicarded_reads', 'SPAdes_number_contigs', 'SPAdes_number_bp', 'SPAdes_filtered_contigs', 'SPAdes_filtered_bp', 'assembly_coverage_initial', 'assembly_coverage_filtered', 'mapped_reads_percentage', 'mapping_filtered_contigs', 'mapping_filtered_bp', 'Pilon_changes', 'Pilon_contigs_changed', 'MLST_scheme', 'MLST_ST', 'final_assembly']

    for directory in directories:
        sample = directory
        directory = os.path.join(os.path.join(innucaOut, directory, ''))

        files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]
        if len(files) == 0:
            print 'No files found! Continue to the next sample'
        else:
            results[sample] = {}
            for field in fields:
                results[sample][field] = 'NA'

            results[sample]['#samples'] = sample

            for file_found in files:
                name_file_found = file_found
                file_found = os.path.join(directory, file_found)

                if name_file_found == 'coverage_report.txt':
                    data = []
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            if len(line) > 0:
                                line = line.splitlines()[0].rsplit('x')[0]
                                data.append(line)
                    if len(data) == 2:
                        results[sample]['first_coverage'] = data[0]
                        results[sample]['second_Coverage'] = data[1]
                    elif len(data) == 1:
                        if samples_report is not None:
                            with open(samples_report, 'rtU') as reader:
                                header = None
                                for line in reader:
                                    line = line.splitlines()[0]
                                    if len(line) > 0:
                                        line = line.split('\t')

                                        if line[0].startswith('#'):
                                            header = line

                                        if line[0] in sample:
                                            if line[header.index("first_Coverage_runSuccessfully")] == 'True':
                                                results[sample]['first_coverage'] = data[0]
                                            elif line[header.index("second_Coverage_runSuccessfully")] == 'True':
                                                results[sample]['second_Coverage'] = data[0]
                elif name_file_found == 'trueCoverage_report.txt':
                    general, absent_genes, multiple_alleles, sample_coverage = False, False, False, False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    if line.startswith('general', 1):
                                        general = True
                                    else:
                                        general = False
                                elif line.startswith('>') and general:
                                    if line.startswith('number_absent_genes', 1):
                                        absent_genes = True
                                    elif line.startswith('number_genes_multiple_alleles', 1):
                                        multiple_alleles = True
                                    elif line.startswith('mean_sample_coverage', 1):
                                        sample_coverage = True
                                else:
                                    if general:
                                        if absent_genes:
                                            results[sample]['trueCoverage_absent_genes'] = line
                                            absent_genes = False
                                        elif multiple_alleles:
                                            results[sample]['trueCoverage_multiple_alleles'] = line
                                            multiple_alleles = False
                                        elif sample_coverage:
                                            results[sample]['trueCoverage_sample_coverage'] = line
                                            sample_coverage = False
                elif name_file_found == 'pear_report.txt':
                    assembled, unassembled, discarded = False, False, False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    if line.startswith('assembled_reads', 1):
                                        assembled = True
                                    elif line.startswith('unassembled_reads', 1):
                                        unassembled = True
                                    elif line.startswith('discarded_reads', 1):
                                        discarded = True
                                else:
                                    if assembled:
                                        results[sample]['pear_assembled_reads'] = line
                                        assembled = False
                                    elif unassembled:
                                        results[sample]['pear_unassembled_reads'] = line
                                        unassembled = False
                                    elif discarded:
                                        results[sample]['pear_discarded_reads'] = line
                                        discarded = False
                elif name_file_found.startswith('spades_report.original.'):
                    general = False

                    contigs = False
                    bp = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            if len(line) > 0:
                                line = line.splitlines()[0].split(' ')[0]
                                if line.startswith('#'):
                                    if line.startswith('general', 1):
                                        general = True
                                    else:
                                        general = False
                                elif line.startswith('>') and general:
                                    if line.startswith('contigs', 1):
                                        contigs = True
                                    elif line.startswith('bp', 1):
                                        bp = True
                                else:
                                    if contigs:
                                        results[sample]['SPAdes_number_contigs'] = line
                                        contigs = False
                                    if bp:
                                        results[sample]['SPAdes_number_bp'] = line
                                        bp = False
                elif name_file_found.startswith('spades_report.filtered.'):
                    general = False

                    contigs = False
                    bp = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            if len(line) > 0:
                                line = line.splitlines()[0].split(' ')[0]
                                if line.startswith('#'):
                                    if line.startswith('general', 1):
                                        general = True
                                    else:
                                        general = False
                                elif line.startswith('>') and general:
                                    if line.startswith('contigs', 1):
                                        contigs = True
                                    elif line.startswith('bp', 1):
                                        bp = True
                                else:
                                    if contigs:
                                        results[sample]['SPAdes_filtered_contigs'] = line
                                        contigs = False
                                    elif bp:
                                        results[sample]['SPAdes_filtered_bp'] = line
                                        bp = False
                elif name_file_found == 'assembly_mapping_report.coverage.txt':
                    general = False

                    initial = False
                    filtered = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    if line.startswith('general', 1):
                                        general = True
                                    else:
                                        general = False
                                elif line.startswith('>') and general:
                                    if line.startswith('initial', 1):
                                        initial = True
                                    elif line.startswith('filtered', 1):
                                        filtered = True
                                else:
                                    if initial:
                                        results[sample]['assembly_coverage_initial'] = line
                                        initial = False
                                    elif filtered:
                                        results[sample]['assembly_coverage_filtered'] = line
                                        filtered = False
                elif name_file_found == 'assembly_mapping_report.mapping.txt':
                    total = False
                    total_reads = 0
                    mapped = False
                    mapped_reads = 0
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0].split(' ')[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    total = False
                                    mapped = False
                                    if line.startswith('in_total', 1):
                                        total = True
                                    elif line.startswith('mapped', 1):
                                        mapped = True
                                elif line.startswith('>'):
                                    continue
                                else:
                                    if total:
                                        total_reads = int(line)
                                        total = False
                                    elif mapped:
                                        mapped_reads = int(line)
                                        mapped = False
                    results[sample]['mapped_reads_percentage'] = str(round((float(mapped_reads) / float(total_reads)), 2) * 100)
                elif name_file_found.startswith('assembly_mapping_report.sequences_filtered.'):
                    general = False

                    contigs = False
                    bp = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    if line.startswith('general', 1):
                                        general = True
                                    else:
                                        general = False
                                elif line.startswith('>') and general:
                                    if line.startswith('contigs', 1):
                                        contigs = True
                                    elif line.startswith('bp', 1):
                                        bp = True
                                else:
                                    if contigs:
                                        results[sample]['mapping_filtered_contigs'] = line
                                        contigs = False
                                    if bp:
                                        results[sample]['mapping_filtered_bp'] = line
                                        bp = False
                elif name_file_found == 'pilon_report.txt':
                    general = False

                    changes = False
                    contigs = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            if len(line) > 0:
                                line = line.splitlines()[0]
                                if line != '#by_contigs':
                                    if line.startswith('#'):
                                        if line.startswith('general', 1):
                                            general = True
                                        else:
                                            general = False
                                    elif line.startswith('>') and general:
                                        if line.startswith('changes', 1):
                                            changes = True
                                        elif line.startswith('contigs', 1):
                                            contigs = True
                                    else:
                                        if changes:
                                            results[sample]['Pilon_changes'] = line
                                            changes = False
                                        if contigs:
                                            results[sample]['Pilon_contigs_changed'] = line
                                            contigs = False
                                else:
                                    break
                elif name_file_found == 'mlst_report.txt':
                    species = False
                    st = False
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                if line.startswith('#'):
                                    if line.startswith('scheme', 1):
                                        species = True
                                    elif line.startswith('ST', 1):
                                        st = True
                                else:
                                    if species:
                                        results[sample]['MLST_scheme'] = line
                                        species = False
                                    if st:
                                        results[sample]['MLST_ST'] = line
                                        st = False
                elif name_file_found == 'final_assembly.txt':
                    with open(file_found, 'rtU') as reader:
                        for line in reader:
                            line = line.splitlines()[0]
                            if len(line) > 0:
                                results[sample]['final_assembly'] = line

    if len(results) == 0:
        sys.exit('No results were found')

    print '\n' + 'Writing results...'
    with open(os.path.join(outdir, str('combine_samples_reports.' + time_str + '.tab')), 'wt') as report:
        report.write('\t'.join(fields) + '\n')
        for sample in results:
            sample_data = [results[sample][field] for field in fields]
            report.write('\t'.join(sample_data) + '\n')

    if json:
        import json
        with open(os.path.join(outdir, str('combine_samples_reports.' + time_str + '.json')), 'wt') as writer:
            json.dump(results, writer)

    if number_samples is None:
        number_samples = len(results)
    print '\n' + 'DONE: {number_samples} samples analysed.'.format(number_samples=number_samples)


def check_create_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


def main():

    parser = argparse.ArgumentParser(prog='python combine_reports.py', description="Combine INNUca reports (Estimated Coverage, True Coverage, Pear, SPAdes, Pilon, Assembly Mapping, MLST)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--innucaOut', type=str, metavar='/path/to/INNUca/output/directory/', help='Path to INNUca output directory', required=True)

    parser_optional = parser.add_argument_group('Facultative options')
    parser_optional.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/', help='Path to where to store the outputs', required=False, default='.')
    parser_optional.add_argument('--json', action='store_true', help='Also save the results in json format')

    args = parser.parse_args()

    combine_reports(args.innucaOut, args.outdir, args.json, None, None)


if __name__ == "__main__":
    main()
