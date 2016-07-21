INNUca.py
=========
INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search*

<https://github.com/B-UMMI/INNUca>

Requirements
------------

 - Illumina paired-end reads
 - Expected species name
 - Expected genome size in Mb

Dependencies
------------
**Mandatory**

 - *Java JDK* >= v1.8
 - *mlst* (https://github.com/tseemann/mlst) >= v2.4 (it is recommended
   to use a mlst version with updated databases)
 - *gzip* >= v1.6 (normally found in Linux OS)
 - *bzip2* >= v1.0.6 (normally found in Linux OS)

**Optional**
(executables are provided, but user's own executables can be used with `--doNotUseProvidedSoftware` option)

 - *FastQC* = v0.11.5
 - *Trimmomatic* = v0.36 (make sure the .jar file is executable and it is
   in your PATH)
 - *SPAdes* >= v3.7.1

Installation
------------
    git clone https://github.com/B-UMMI/INNUca.git

Usage
-----

    usage: INNUca.py [-h] [--version]
                     -i /path/to/input/directory/
                     -s "Streptococcus agalactiae"
                     -g 2.1
                     [-o /output/directory/] [-j N] [--doNotUseProvidedSoftware]
                     [--pairEnd_filesSeparation "_left/rigth.fq.gz" "_left/rigth.fq.gz"]
                     [--skipEstimatedCoverage] [--skipFastQC] [--skipTrimmomatic]
                     [--skipSPAdes] [--skipMLST]
                     [--adapters adaptersFile.fasta | --doNotSearchAdapters]
                     [--doNotTrimCrops | [[--trimCrop N] [--trimHeadCrop N]]]
                     [--trimSlidingWindow window:meanQuality] [--trimMinLength N]
                     [--trimLeading N] [--trimTrailing N]
                     [--spadesNotUseCareful] [--spadesMinContigsLength N]
                     [--spadesKmers 55,77] [--spadesMaxMemory N]
                     [--spadesMinCoverage 10]

    INNUca - Reads Control and Assembly

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    Required options:
      -i /path/to/input/directory/, --inputDirectory /path/to/input/directory/
                            Path to directory containing the fastq files. Can be
                            organized in separete directories by samples or all
                            together (default: None)
      -s "Streptococcus agalactiae", --speciesExpected "Streptococcus agalactiae"
                            Expected species name (default: None)
      -g 2.1, --genomeSizeExpectedMb 2.1
                            Expected genome size in Mb (default: None)

    General options:
      -o /output/directory/, --outdir /output/directory/
                            Path for output directory (default: .)
      -j N, --threads N     Number of threads (default: [1])
      --doNotUseProvidedSoftware
                            Tells the software to not use FastQC, Trimmomatic,
                            SPAdes and Samtools that are provided with INNUca.py
                            (default: False)
      --pairEnd_filesSeparation "_left/rigth.fq.gz" "_left/rigth.fq.gz"
                            For unusual pair-end files separation designations,
                            you can provide two strings containning the end of
                            fastq files names to designate each file from a pair-
                            end data ("_left.fq.gz" "_rigth.fq.gz" for
                            sample_left.fq.gz sample_right.fq.gz) (default: None)
      --skipEstimatedCoverage
                            Tells the programme to not estimate coverage depth
                            based on number of sequenced nucleotides and expected
                            genome size (default: False)
      --skipFastQC          Tells the programme to not run FastQC analysis
                            (default: False)
      --skipTrimmomatic     Tells the programme to not run Trimmomatic (default:
                            False)
      --skipSPAdes          Tells the programme to not run SPAdes and consequently
                            MLST analysis (requires SPAdes contigs) (default:
                            False)
      --skipMLST            Tells the programme to not run MLST analysis (default:
                            False)

    Adapters options (one of the following):
      --adapters adaptersFile.fasta
                            Fasta file containing adapters sequences to be used in
                            FastQC and Trimmomatic (default: [None])
      --doNotSearchAdapters
                            Tells INNUca.py to not search for adapters and clip
                            them during Trimmomatic step (default: False)

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
                            falls below a threshold (default: ['5:20'])
      --trimLeading N       Trimmomatic: cut bases off the start of a read, if
                            below a threshold quality (default: [3])
      --trimTrailing N      Trimmomatic: cut bases off the end of a read, if below
                            a threshold quality (default: [3])
      --trimMinLength N     Trimmomatic: drop the read if it is below a specified
                            length (default: [55])
      --trimKeepFiles       Tells INNUca.py to not remove the output of Trimmomatic

    SPAdes options:
      --spadesNotUseCareful
                            Tells SPAdes to only perform the assembly without the
                            --careful option (default: False)
      --spadesMinContigsLength N
                            Filter SPAdes contigs for length greater or equal than
                            this value (default: [200])
      --spadesKmers 55,77   Manually sets SPAdes k-mers lengths (all values must
                            be odd, less than 128) (default: [55, 77, 99, 113,
                            127])
      --spadesMaxMemory N   The maximum amount of RAM Gb for SPAdes to use
                            (default: [25])
      --spadesMinCoverage 10
                            The minimum number of reads to consider an edge in the  
                            de Bruijn graph (or path I am not sure). Can also be
                            auto or off (default: ['off'])
      --spadesSaveReport    Tells INNUca to store the number of contigs and
                            assembled nucleotides for each sample

Contact
-------
Miguel Machado
<mpmachado@medicina.ulisboa.pt>



> Written with [StackEdit](https://stackedit.io/).
