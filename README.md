# INNUca.py

INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality
assessment, and possible contamination detection*

<https://github.com/B-UMMI/INNUca>

---

## Requirements

 - Illumina paired-end reads (paired end information: sampleName_R1_001 / sampleName_R2_001 OR sampleName_1 / sampleName_2) (gzip compressed: .fastq.gz or .fq.gz)
 - Expected species name
 - Expected genome size in Mb

## Dependencies

**Mandatory**

 - *Java JDK*
 - [*Kraken*](https://ccb.jhu.edu/software/kraken/) >= v0.10.6 with *Kraken* DB (whenever *Kraken* module should run)
 - [*mlst*](https://github.com/tseemann/mlst) >= v2.4 (whenever *mlst* module should run) (it is recommended to use a mlst version with updated databases)
 - [*ReMatCh*](https://github.com/B-UMMI/ReMatCh) >= v3.2 (whenever *true coverage* module should run)
 - *gzip* >= v1.6 (normally found in Linux OS)

**Optional**
(executables are provided, but user's own executables can be used with `--doNotUseProvidedSoftware` option)

 - *Bowtie2* >= v2.2.9
 - *Samtools* = v1.3.1
 - *FastQC* = v0.11.5
 - *Trimmomatic* = v0.36 (make sure the .jar file is executable and it is
   in your PATH)
 - *Pear* = v0.9.10
 - *SPAdes* >= v3.9.0
 - *Pilon* = v1.18

## Installation

    git clone https://github.com/B-UMMI/INNUca.git

---

## Usage

```
usage: INNUca.py [-h] [--version] -s "Streptococcus agalactiae" -g 2.1
                 (-i /path/to/input/directory/ | -f /path/to/input/file_1.fq.gz /path/to/input/file_2.fq.gz)
                 [-o /output/directory/] [-j N] [--keepBAM] [--noLog]
                 [--noGitInfo] [--json] [--jarMaxMemory 10]
                 [--doNotUseProvidedSoftware]
                 [--runKraken] [--skipEstimatedCoverage] [--skipTrueCoverage]
                 [--skipFastQC] [--skipTrimmomatic] [--runPear] [--skipSPAdes]
                 [--skipAssemblyMapping] [--skipPilon] [--skipMLST]
                 [--adapters adaptersFile.fasta | --doNotSearchAdapters]
                 [--krakenDB minikraken_20171013_4GB] [--krakenQuick]
                 [--krakenProceed] [--krakenIgnoreQC] [--krakenMemory]
                 [--krakenMinCov 1.5] [--krakenMaxUnclass 1.5]
                 [--krakenMinQual N]
                 [--estimatedMinimumCoverage N]
                 [--trueConfigFile species.config]
                 [--fastQCkeepFiles] [--fastQCproceed]
                 [--trimKeepFiles] [--doNotTrimCrops] [--trimCrop N]
                 [--trimHeadCrop N] [--trimSlidingWindow window:meanQuality]
                 [--trimLeading N] [--trimTrailing N] [--trimMinLength N]
                 [--spadesVersion 3.11.0] [--spadesNotUseCareful]
                 [--spadesMinContigsLength N] [--spadesMaxMemory N]
                 [--spadesMinCoverageAssembly N] [--spadesMinKmerCovContigs N]
                 [--spadesKmers 55 77 [55 77 ...] | --spadesDefaultKmers]
                 [--assemblyMinCoverageContigs N] [--maxNumberContigs N]
                 [--saveExcludedContigs] [--keepIntermediateAssemblies]
                 [--pilonKeepFiles]
                 [--mlstIgnoreQC]
                 [--pearKeepFiles] [--pearMinOverlap N]

INNUca - Reads Control and Assembly

optional arguments:
  -h, --help            show this help message and exit
  --version             Version information

Required options:
  -s "Streptococcus agalactiae", --speciesExpected "Streptococcus agalactiae"
                        Expected species name (default: None)
  -g 2.1, --genomeSizeExpectedMb 2.1
                        Expected genome size in Mb (default: None)

Required INPUT options (one of the following):
  -i /path/to/input/directory/, --inputDirectory /path/to/input/directory/
                        Path to directory containing the fastq files. Can be
                        organized in separete directories by samples or all
                        together (default: None)
  -f /path/to/input/file_1.fq.gz /path/to/input/file_2.fq.gz, --fastq /path/to/input/file_1.fq.gz /path/to/input/file_2.fq.gz
                        Path to Pair-End Fastq files (default: None)

General options:
  -o /output/directory/, --outdir /output/directory/
                        Path for output directory (default: ./) (default: .)
  -j N, --threads N     Number of threads (default: 1) (default: 1)
  --keepBAM             Keep the last BAM file produced (with mapped and
                        unmapped reads) (default: False)
  --noLog               Do not create a log file (default: False)
  --noGitInfo           Do not retrieve GitHub repository information
                        (default: False)
  --json                Tells INNUca to save the results also in json format
                        (default: False)
  --jarMaxMemory 10     Sets the maximum RAM Gb usage by jar files
                        (Trimmomatic and Pilon). Can also be auto or off. When
                        auto is set, 1 Gb per thread will be used up to the
                        free available memory. (default: off (default: off)
  --doNotUseProvidedSoftware
                        Tells the software to not use FastQC, Trimmomatic,
                        SPAdes, Bowtie2, Samtools and Pilon that are provided
                        with INNUca.py (default: False)

Running modules options:
  --runKraken           Sets INNUca to run Kraken (default: False)
  --skipEstimatedCoverage
                        Tells the programme to not estimate coverage depth
                        based on number of sequenced nucleotides and expected
                        genome size (default: False)
  --skipTrueCoverage    Tells the programme to not run trueCoverage_ReMatCh
                        analysis (default: False)
  --skipFastQC          Tells the programme to not run FastQC analysis
                        (default: False)
  --skipTrimmomatic     Tells the programme to not run Trimmomatic (default:
                        False)
  --runPear             Tells the programme to run Pear (default: False)
  --skipSPAdes          Tells the programme to not run SPAdes and consequently
                        all the modules that are assembly based (Assembly
                        Mapping check, Pilon correction, MLST analysis and
                        Kraken on the assembly) (default: False)
  --skipAssemblyMapping
                        Tells the programme to not run Assembly Mapping check
                        (default: False)
  --skipPilon           Tells the programme to not run Pilon correction
                        (default: False)
  --skipMLST            Tells the programme to not run MLST analysis (default:
                        False)

Adapters options (one of the following):
  Control how adapters are handle by INNUca. If none of these options are provided, INNUca
  will use the Nextera XT and PE TruSeq files found at INNUca/src/Trimmomatic-0.36/adapters/

  --adapters adaptersFile.fasta
                        Fasta file containing adapters sequences to be used in
                        FastQC and Trimmomatic (default: None)
  --doNotSearchAdapters
                        Tells INNUca.py to not search for adapters and clip
                        them during Trimmomatic step (default: False)

Kraken options:
  --krakenDB minikraken_20171013_4GB
                        Name of Kraken DB found in path, or complete path to
                        the directory containing the Kraken DB files (for
                        example /path/to/directory/minikraken_20171013_4GB)
                        (default: None)
  --krakenQuick         Set Kraken to do a quick operation and only use the
                        first hits (default: False)
  --krakenProceed       Do not stop INNUca.py if sample fails Kraken (default:
                        False)
  --krakenIgnoreQC      Ignore Kraken QA/QC in sample quality assessment.
                        Useful when analysing data from possible new species
                        or higher taxonomic levels (higher than species)
                        (default: False)
  --krakenMemory        Set Kraken to load the DB into the memory before run
                        (default: False)
  --krakenMinCov 1.5    Minimum percentage of fragments covered to consider
                        the taxon. If nothing is specified, the hundredth of
                        the taxon found (species or genus if no species are
                        available for a given genus) with higher percentage of
                        fragments covered (excluding unclassified category)
                        will be used. (default: None)
  --krakenMaxUnclass 1.5
                        Maximum percentage of unclassified fragments allowed.
                        If nothing is specified, the tenth of 100 minus the
                        percentage of fragments of the taxon found (species or
                        genus if no species are available for a given genus)
                        with higher percentage of fragments covered (excluding
                        unclassified category) will be used. (default: None)
  --krakenMinQual N     Sets the minimum base quality to be used in
                        classification (default: 10) (default: 10)

Estimated Coverage options:
  This module estimates the depth of coverage by dividing the number of
  sequenced nucleotides (raw or processed reads) by the expected genome size
  (in bps)

  --estimatedMinimumCoverage N
                        Minimum estimated coverage to continue INNUca pipeline
                        (default: 15) (default: 15)

trueCoverage_ReMatCh options:
  This module calculates an improved estimation of the true bacterial
  chromosome coverage via read mapping against reference gene sequences
  distributed throughout the genome. This approach alleviates coverage
  estimation bias introduced by mobile genetic elements and other similar
  occurrences. Moreover, this module can also detect multiple strains or
  species contamination by searching for heterozygous positions. INNUca
  provides target sequences for some species together with the desired
  settings and the QA/QC decision rules (in
  INNUca/modules/trueCoverage_rematch/). If the expected species matches any
  of the species provided files, trueCoverage_ReMatCh module will run.

  --trueConfigFile species.config
                        File with trueCoverage_ReMatCh settings. Some species
                        specific config files can be found in
                        INNUca/modules/trueCoverage_rematch/ folder. Use those
                        files as example files. For species with config files
                        in INNUca/modules/trueCoverage_rematch/ folder (not
                        pre releases versions, marked with "pre."),
                        trueCoverage_ReMatCh will run by default, unless
                        --skipTrueCoverage is specified. Do not use together
                        with --skipTrueCoverage option (default: None)

FastQC options:
  --fastQCkeepFiles     Tells INNUca.py to not remove the output of FastQC
                        (default: False)
  --fastQCproceed       Do not stop INNUca.py if sample fails FastQC (default:
                        False)

Trimmomatic options:
  --trimKeepFiles       Tells INNUca.py to not remove the output of
                        Trimmomatic (default: False)
  --doNotTrimCrops      Tells INNUca.py to not cut the beginning and end of
                        reads during Trimmomatic step (unless specified with
                        --trimCrop or --trimHeadCrop, INNUca.py will search
                        for nucleotide content bias at both ends and will cut
                        by there) (default: False)
  --trimCrop N          Cut the specified number of bases to the end of the
                        maximum reads length. By default, the number of bases
                        to cut is calculated using sample FastQC results and
                        based on the G/C content. The first position of the
                        second half of the reads with GC bias (80 < GC
                        percentage > 120) followed by at least two other
                        biased positions, are marked to be trimmed. (default:
                        None)
  --trimHeadCrop N      Trimmomatic: cut the specified number of bases from
                        the start of the reads. By default, the number of
                        bases to cut is calculated using sample FastQC results
                        and based on the G/C content. The first position of
                        the first half of the reads with GC bias (80 < GC
                        percentage > 120) followed by at least two other
                        unbiased positions, are marked to be trimmed.
                        (default: None)
  --trimSlidingWindow window:meanQuality
                        Trimmomatic: perform a sliding window trimming,
                        cutting once the average quality within the window
                        falls below a threshold (default: 5:20) (default:
                        5:20)
  --trimLeading N       Trimmomatic: cut bases off the start of a read, if
                        below a threshold quality (default: 3) (default: 3)
  --trimTrailing N      Trimmomatic: cut bases off the end of a read, if below
                        a threshold quality (default: 3) (default: 3)
  --trimMinLength N     Trimmomatic: drop the read if it is below a specified
                        length (default: 55) (default: 55)

SPAdes options:
  --spadesVersion 3.11.0
                        Tells INNUca.py which SPAdes version to use (available
                        options: 3.9.0, 3.10.1, 3.11.0) (default: 3.11.0)
                        (default: 3.11.0)
  --spadesNotUseCareful
                        Tells SPAdes to only perform the assembly without the
                        --careful option (default: False)
  --spadesMinContigsLength N
                        Filter SPAdes contigs for length greater or equal than
                        this value (default: maximum reads size or 200 bp)
                        (default: None)
  --spadesMaxMemory N   The maximum amount of RAM Gb for SPAdes to use
                        (default: 2 Gb per thread will be used up to the free
                        available memory) (default: None)
  --spadesMinCoverageAssembly N
                        The minimum number of reads to consider an edge in the
                        de Bruijn graph during the assembly. Can also be auto
                        or off (default: 2) (default: 2)
  --spadesMinKmerCovContigs N
                        Minimum contigs K-mer coverage. After assembly only
                        keep contigs with reported k-mer coverage equal or
                        above this value (default: 2) (default: 2)

SPAdes k-mers options (one of the following):
  --spadesKmers 55 77 [55 77 ...]
                        Manually sets SPAdes k-mers lengths (all values must
                        be odd, lower than 128) (default values: reads length
                        >= 175 [55, 77, 99, 113, 127]; reads length < 175 [21,
                        33, 55, 67, 77]) (default: None)
  --spadesDefaultKmers  Tells INNUca to use SPAdes default k-mers (default:
                        False)

Assembly Mapping options:
  --assemblyMinCoverageContigs N
                        Minimum contigs average coverage. After mapping reads
                        back to the contigs, only keep contigs with at least
                        this average coverage (default: 1/3 of the assembly
                        mean coverage or 10x) (default: None)

Assembly options:
  --maxNumberContigs N  Maximum number of contigs per 1.5 Mb of expected
                        genome size (default: 100) (default: 100)
  --saveExcludedContigs
                        Tells INNUca.py to save excluded contigs (default:
                        False)
  --keepIntermediateAssemblies
                        Tells INNUca to keep all the intermediate assemblies
                        (default: False)

Pilon options:
  --pilonKeepFiles      Tells INNUca.py to not remove the output of Pilon
                        (default: False)

MLST options:
  --mlstIgnoreQC        Ignore MLST QA/QC in sample quality assessment. Useful
                        when analysing data from possible new species or
                        higher taxonomic levels (higher than species)
                        (default: False)

Pear options:
  --pearKeepFiles       Tells INNUca.py to not remove the output of Pear
                        (default: False)
  --pearMinOverlap N    Minimum nucleotide overlap between read pairs for Pear
                        assembly them into only one read (default: 2/3 of
                        maximum reads length determine using FastQC, or
                        Trimmomatic minimum reads length if it runs, or 33
                        nts) (default: None)
```

---


![bioinformaticsopendays_braga_2017_innuca_poster_25](https://user-images.githubusercontent.com/13034956/33260026-9790da3a-d356-11e7-92f0-945dcd0bcec1.png)  
<div style="text-align: right">INNUca's poster presented at **_Bioinformatics Open Days 2017_**, Braga, Portugal (February 23-24)</div>

---

Combine INNUca reports
----------------------
In order to combine **INNUca** reports (Estimate Coverage, True Coverage, Pear, SPAdes, Assembly Mapping, Pilon, MLST), use *combine_reports.py* found in **INNUca** modules folder

    usage: python combine_reports.py [-h] [--version] -i
                              /path/to/INNUca/output/directory/
                              [-o /path/to/output/directory/]

    Combine INNUca reports (Estimated Coverage, True Coverage, Pear, SPAdes, Assembly
    Mapping, Pilon, MLST)

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    Required options:
      -i /path/to/INNUca/output/directory/, --innucaOut /path/to/INNUca/output/directory/
                            Path to INNUca output directory (default: None)

    Facultative options:
      -o /path/to/output/directory/, --outdir /path/to/output/directory/
                            Path to where to store the outputs (default: ['.'])



Combine trueCoverage_ReMatCh module reports
----------------------
In order to manually combine **INNUca** trueCoverage_ReMatCh module reports in respect to gene information, use *combine_trueCoverage_reports.py* found in **INNUca** modules/trueCoverage_rematch folder

    usage: python combine_trueCoverage_reports.py [-h] [--version] -i
                                                  /path/to/INNUca/output/directory/
                                                  [-o /path/to/output/directory/]
                                                  [--minimum_gene_coverage 80]

    Combine trueCoverage_ReMatCh module reports in respect to gene information.

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    Required options:
      -i /path/to/INNUca/output/directory/, --innucaOut /path/to/INNUca/output/directory/
                            Path to INNUca output directory (default: None)

    Facultative options:
      -o /path/to/output/directory/, --outdir /path/to/output/directory/
                            Path to where to store the outputs (default: .)
      --minimum_gene_coverage 80
                            Minimum percentage of sequence length (with a minimum
                            of read depth to consider a position to be present) to
                            determine whether a gene is present. (default: 80)



## Citation
MP Machado, J Halkilahti, A Jaakkonen, DN Silva, I Mendes, Y Nalbantoglu, V Borges, M Ramirez, M Rossi, JA Carriço. _INNUca_ **GitHub** https://github.com/B-UMMI/INNUca

Contact
-------
Miguel Machado
<mpmachado@medicina.ulisboa.pt>



> Written with [StackEdit](https://stackedit.io/).
