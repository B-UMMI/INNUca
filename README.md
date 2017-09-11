INNUca.py
=========
INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search*

<https://github.com/B-UMMI/INNUca>

Requirements
------------

 - Illumina paired-end reads (paired end information: sampleName_R1_001 / sampleName_R2_001 OR sampleName_1 / sampleName_2) (gzip compressed: .fastq.gz or .fq.gz)
 - Expected species name
 - Expected genome size in Mb

Dependencies
------------
**Mandatory**

 - *Java JDK*
 - *mlst* (https://github.com/tseemann/mlst) >= v2.4 (whenever *mlst* module should run) (it is recommended to use a mlst version with updated databases)
 - *ReMatCh* (https://github.com/B-UMMI/ReMatCh) >= v3.2 (whenever *true coverage* module should run)
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

Installation
------------
    git clone https://github.com/B-UMMI/INNUca.git

Usage
-----

    usage: INNUca.py [-h] [--version]
                     -s "Streptococcus agalactiae" -g 2.1
                     (-i /path/to/input/directory/ | -f /path/to/input/file_1.fq.gz /path/to/input/file_2.fq.gz)
                     [-o /output/directory/] [-j N]
                     [--jarMaxMemory 10] [--doNotUseProvidedSoftware]
                     [--keepIntermediateAssemblies]
                     [--skipEstimatedCoverage] [--skipFastQC]
                     [--skipTrimmomatic] [--skipSPAdes] [--skipAssemblyMapping]
                     [--skipPilon] [--skipMLST] [--runPear] [--noLog] [--noGitInfo]
                     [--skipTrueCoverage | --trueConfigFile species.config]
                     [--adapters adaptersFile.fasta | --doNotSearchAdapters]
                     [--estimatedMinimumCoverage N]
                     [--fastQCkeepFiles] [--fastQCproceed]
                     [--doNotTrimCrops | [[--trimCrop N] [--trimHeadCrop N]]]
                     [--trimSlidingWindow window:meanQuality] [--trimLeading N]
                     [--trimTrailing N] [--trimMinLength N] [--trimKeepFiles]
                     [--pearKeepFiles] [--pearMinOverlap N]
                     [--spadesVersion] [--spadesNotUseCareful]
                     [--spadesMinContigsLength N] [--spadesMaxMemory N]
                     [--spadesMinCoverageAssembly 10] [--spadesMinKmerCovContigs N]
                     [--spadesKmers 55 77 [55 77 ...] | --spadesDefaultKmers]
                     [--assemblyMinCoverageContigs N]
                     [--maxNumberContigs N] [--saveExcludedContigs]
                     [--pilonKeepFiles]

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
                            Path for output directory (default: .)
      -j N, --threads N     Number of threads (default: 1)
      --jarMaxMemory 10     Sets the maximum RAM Gb usage by jar files
                            (Trimmomatic and Pilon). Can also be auto or off. When
                            auto is set, 1 Gb per thread will be used up to the
                            free available memory (default: off)
      --doNotUseProvidedSoftware
                            Tells the software to not use FastQC, Trimmomatic,
                            SPAdes, Bowtie2, Samtools and Pilon that are provided
                            with INNUca.py (default: False)
      --keepIntermediateAssemblies
                            Tells INNUca to keep all the intermediate assemblies
                            (default: False)
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
      --skipSPAdes          Tells the programme to not run SPAdes and consequently
                            Pilon correction, Assembly Mapping check and MLST
                            analysis (SPAdes contigs required) (default: False)
      --skipAssemblyMapping
                            Tells the programme to not run Assembly Mapping check
                            (default: False)
      --skipPilon           Tells the programme to not run Pilon correction and
                            consequently Assembly Mapping check (bam files
                            required) (default: False)
      --skipMLST            Tells the programme to not run MLST analysis (default:
                            False)
      --runPear             Tells the programme to run Pear (default: False)
      --noLog               Do not create a log file (default: False)
      --noGitInfo           Do not retreive GitHub repository information
                            (default: False)

    Adapters options (one of the following):
      --adapters adaptersFile.fasta
                            Fasta file containing adapters sequences to be used in
                            FastQC and Trimmomatic (default: None)
      --doNotSearchAdapters
                            Tells INNUca.py to not search for adapters and clip
                            them during Trimmomatic step (default: False)

    Estimated Coverage options:
      --estimatedMinimumCoverage N
                            Minimum estimated coverage to continue INNUca pipeline
                            (default: 15)

    trueCoverage_ReMatCh options:
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
      --fastQCkeepFiles     Tells INNUca.py to not remove the output of
                            FastQC (default: False)
      --fastQCproceed       Do not stop INNUca.py if sample fails FastQC (default:
                            False)

    Trimmomatic options:
      --doNotTrimCrops      Tells INNUca.py to not cut the beginning and end of
                            reads during Trimmomatic step (unless specified with
                            --trimCrop or --trimHeadCrop, INNUca.py will search
                            for nucleotide content bias at both ends and will cut
                            by there) (default: False)
      --trimCrop N          Cut the specified number of bases to the end of the
                            maximum reads length (default: None)
      --trimHeadCrop N      Trimmomatic: cut the specified number of bases from
                            the start of the reads (default: None)
      --trimSlidingWindow window:meanQuality
                            Trimmomatic: perform a sliding window trimming,
                            cutting once the average quality within the window
                            falls below a threshold (default: 5:20)
      --trimLeading N       Trimmomatic: cut bases off the start of a read, if
                            below a threshold quality (default: 3)
      --trimTrailing N      Trimmomatic: cut bases off the end of a read, if below
                            a threshold quality (default: 3)
      --trimMinLength N     Trimmomatic: drop the read if it is below a specified
                            length (default: 55)
      --trimKeepFiles       Tells INNUca.py to not remove the output of
                            Trimmomatic (default: False)

    Pear options:
      --pearKeepFiles       Tells INNUca.py to not remove the output of Pear
                            (default: False)
      --pearMinOverlap      Minimum nucleotide overlap between read pairs for Pear
                            assembly them into only one read (default: 2/3 of maximum
                            reads length or 33 whenever is was not possible to determine
                            it with FastQC)

    SPAdes options:
      --spadesVersion 3.11.0
                            Tells INNUca.py which SPAdes version to use
                            (available options: 3.9.0, 3.10.1, 3.11.0) (default:
                            3.11.0)
      --spadesNotUseCareful
                            Tells SPAdes to only perform the assembly without the
                            --careful option (default: False)
      --spadesMinContigsLength N
                            Filter SPAdes contigs for length greater or equal than
                            this value (default: maximum reads size or 200 bp)
      --spadesMaxMemory N   The maximum amount of RAM Gb for SPAdes to use
                            (default: 2 Gb per thread will be used up to the free
                            available memory)
      --spadesMinCoverageAssembly 10
                            The minimum number of reads to consider an edge in the
                            de Bruijn graph during the assembly. Can also be auto
                            or off (default: 2)
      --spadesMinKmerCovContigs N
                            Minimum contigs K-mer coverage. After assembly only
                            keep contigs with reported k-mer coverage equal or
                            above this value (default: 2)

    SPAdes k-mers options (one of the following):
      --spadesKmers 55 77 [55 77 ...]
                            Manually sets SPAdes k-mers lengths (all values must
                            be odd, lower than 128) (default values: reads
                            length >= 175 [55, 77, 99, 113, 127]; reads
                            length < 175 [21, 33, 55, 67, 77])
      --spadesDefaultKmers  Tells INNUca to use SPAdes default k-mers (default:
                            False)

    Assembly Mapping options:
      --assemblyMinCoverageContigs N
                            Minimum contigs average coverage. After mapping reads
                            back to the contigs, only keep contigs with at least
                            this average coverage (default: 1/3 of the assembly
                            mean coverage or 10x)

    Assembly options:
      --maxNumberContigs N  Maximum number of contigs per 1.5 Mb of expected
                            genome size (default: 100)
      --saveExcludedContigs Tells INNUca.py to save excluded contigs (default: False)

    Pilon options:
      --pilonKeepFiles      Tells INNUca.py to not remove the output of Pilon
                            (default: False)



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



Contact
-------
Miguel Machado
<mpmachado@medicina.ulisboa.pt>



> Written with [StackEdit](https://stackedit.io/).
