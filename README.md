### See `docs/index.md` for our wiki documentation.

# Team1-GenomeAssembly
Mariam Nawaz, Pranoti Harkud, Logan Gloster, Zun Wang, Zack Mudge, Mannan Bhola, Yi-Ming Chen, Palak Aggarwal 

## Overview of Assembly Pipeline
The assembly pipeline performs read trimming and de novo assembly given paired-end read files in FASTQ format. Initial read quality assessment and trimming are performed with FastQC and BBDuk, respectively. MultiQC is used to compile read quality metrics for all isolates into one quality assessment output. Assemblies are produced using SPAdes, and assembly quality reports are produced using QUAST.

### Trimming
Trims the left and right ends of each read using BBDuk, with a quality score threshold. The default setting for the assembly script is to start at Q5, and increase by 1 until the reads pass FastQC.
### De novo assembly
SPAdes produces assemblies for different k-mer sizes: 21, 33, 55, 77, 99, and 127. It chooses the optimal assembly from the individual assemblies.
### Assembly Quality
QUAST is used to report assembly quality metrics such as total length, N50, and number of contigs.

## Installing Dependencies
1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Copy `environment.yml` to a local directory. You can change the first line of `environment.yml` to specify the name for your environment.
3. Install dependencies: `conda env create -f environment.yml`
4. `conda activate <name_of_env>`

## Running Assembly Pipeline
### Usage and Options
See "Directory Structure" section below for more details on file placement.
    
       
    usage: assembly.py [-h] [--threshold THRESHOLD] [--skip_trimming] [--skip_quast] data_dir

    positional arguments:
    data_dir            Directory containing the original
                        <isolate_name>_1.fq.gz and <isolate_name>_2.fq.gz
                        files for each isolate OR, if the --skip-trimming
                        option is given, the directory containing the trimmed
                        and unzipped .fq files for each isolate

    optional arguments:
    -h, --help          show this help message and exit
    --threshold THRESHOLD
                        Trimming quality threshold to start at. E.g. 5 starts
                        trimming the read ends at Q5, and increases quality
                        threshold by 1 until trimmed reads pass FASTQC
    --skip_trimming     Skip trimming. If using this option, provide the
                        directory containing the trimmed reads as <data_dir>
    --skip_quast        Skip QUAST of final assemblies.

### Directory Structure    
    . 
    ├── assembly.py
    └───data_dir
        │
        └───CGT1005
        |        CGT1005_1.fq.gz
        |        CGT1005_2.fq.gz
        └───CGT1029
        |        CGT1029_1.fq.gz
        |        CGT1029_2.fq.gz
        ...

### Example: Running assembly pipeline in the background, redirecting stdout to `log.out`, and stderr to `log.err`
`nohup python -u assembly.py data/ > log.out 2> log.err &`

### Output locations
- `assemblies/`: Final assemblies (`contigs.fasta`). Also contains the QUAST outputs for all isolates if QUAST was not skipped.
- `spades_outputs/`: All of the outputs from SPAdes, including assemblies for different k-mer sizes. SPAdes finds the optimal k-mer size based on multiple runs. The assembly for the optimal k-mer size is copied to the `assemblies/` directory at the end.
- `trimmed/`: If trimming was not skipped, FASTQ files with trimmed reads.
- `fastqc/`: If trimming was not skipped, the FastQC outputs for each isolate.
- `multiqc_data/`: If trimming was not skipped, the MultiQC outputs, not including the html report.
- `multiqc_report.html`: If trimming was not skipped, the MultiQC html report.
