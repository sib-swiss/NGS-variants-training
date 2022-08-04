
## Learning outcomes

**After having completed this chapter you will be able to:**

- Describe the general workflow of library preparation and sequencing with an Illumina sequencer
- Explain how the fastq format stores sequence and base quality information
- Calculate probability from phred quality and the other way around
- Explain why base quality and mapping quality are important for detecting variants
- Illustrate the difference between short-read and long-read sequencing
- Explain which type of invention led to development of long-read sequencing
- Explain what impact long read sequencing can have on variant analysis
- Describe how alignment information is stored in a sequence alignment (`.sam`) file
- Define a duplicate alignment and explain how alignment duplicates can affect variant analysis
- Perform an alignment of genomic reads with `bwa mem`
- Generate and interpret the alignment statistics from `samtools flagstat`

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_alignment.pdf){: .md-button }

## Exercises

### 1. Download data and prepare the reference genome

Let's start with the first script of our 'pipeline'. We will use it to download and unpack the course data. Use the code snippet below to create a script called `01_download_course_data.sh`. Store it in `~/workdir/scripts/`, and run it.

```sh title="01_download_course_data.sh"
#!/usr/bin/env bash
cd ~/workdir

wget https://ngs-variants-training.s3.eu-central-1.amazonaws.com/ngs-variants-training.tar.gz
tar -xvf ngs-variants-training.tar.gz
rm ngs-variants-training.tar.gz
```

**Exercise:** This will create the directory `data`. Check out what's in there.

??? done "Answer"
    The directory data contains the following:
    ```
    data
    ├── fastq
    │   ├── father_R1.fastq.gz
    │   ├── father_R2.fastq.gz
    │   ├── mother_R1.fastq.gz
    │   ├── mother_R2.fastq.gz
    │   ├── son_R1.fastq.gz
    │   └── son_R2.fastq.gz
    ├── reference
    │   └── Homo_sapiens.GRCh38.dna.chromosome.20.fa
    └── variants
        ├── 1000g_gold_standard.indels.filtered.vcf
        ├── GCF.38.filtered.renamed.vcf
        ├── NA12878.vcf.gz
        └── NA12878.vcf.gz.tbi

    3 directories, 11 files
    ```

    These are:
    
    * input reads (at `fastq`)
    * a part of the human reference genome (at `reference`)
    * some vcfs with variants for calibration and evaluation (at `variants`)

!!! note "Use `data` only for input"
    The directory `data` that you have just downloaded, contains only input files for the exercises. So, don't write output (except for indexes) to this directory.

In order to index the reference sequence we are going to need some bioinformatics tools. All required tools are pre-installed in a conda environment called `ngs-tools`. In order to use them, every time you open a new terminal, you will have to load the environment:

```sh
conda activate ngs-tools
```

!!! warning "activate ngs-tools"
    Every time you open a new terminal you will have to activate the environment again.

The software `bwa` is in this environment. We will use it for the alignment. Like all alignment software, it requires an index of the reference genome. You can make an index like this:

```sh
bwa index <reference.fa>
```

Make an index of the reference sequence of chromosome 20 of the human genome. You can find the fasta file in `~/workdir/data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa`. Do it with a script called `02_create_bwa_index.sh`. 

??? done "Answer"
    ```sh title="02_create_bwa_index.sh"
    #!/usr/bin/env bash

    cd ~/workdir/data/reference/
    bwa index Homo_sapiens.GRCh38.dna.chromosome.20.fa
    ```

### 2. Read alignment

Check out the [synopsis and manual of `bwa mem`](http://bio-bwa.sourceforge.net/bwa.shtml). We'll be using paired-end reads of three samples that can be found at `~/workdir/data/fastq`. If we run `bwa mem` with default options, which three arguments do we need?

??? done "Answer"
    The manual says:
    ```
    bwa mem [-aCHMpP] [-t nThreads] [-k minSeedLen] ... db.prefix reads.fq [mates.fq]
    ```
    So, we'll need:

    * a database prefix (`db.prefix`)
    * forward reads (`reads.fq`)
    * reverse reads (`mates.fq`)

    For our reference sequence a command would look like:

    ```sh
    cd ~/workdir/

    bwa mem \
    data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    <forward_reads.fq> \
    <reverse_reads.fq> \
    > <alignment.sam>
    ```

We will now go through all the steps concerning alignment for the sample `mother`. To store the results of these steps, we will create a directory within `~/workdir` called `results`. For the alignment, make a script called `03_alignment.sh`. Since we will perform a similar analysis later on for all samples, we store this script in a subdirectory of `~/workdir/scripts` called `mother_only`. 

Long story short, organize your directory `~/workdir/scripts` like this:

```
scripts
├── 01_download_course_data.sh
├── 02_create_bwa_index.sh
└── mother_only
    └── 03_alignment.sh

1 directory, 3 files
```

In `03_alignment.sh` write the commands to perform an alignment with `bwa mem` of the reads from the mother (`mother_R1.fastq` and `mother_R2.fastq`) against chromosome 20. Write the resulting `.sam` file to a directory in `~/workdir/results` called `alignments`.

!!! note "Index prefix is the same a reference filename"
    With default values, the name of the index of a reference for `bwa mem` is the same as the name of the reference itself. In this case, this would be `Homo_sapiens.GRCh38.dna.chromosome.20.fa`.

??? done "Answer"
    Put the script in `~/workdir/scripts/mother_only/`. 

    ```sh title="03_read_alignment.sh"
    #!/usr/bin/env bash

    cd ~/workdir/
    mkdir -p results/alignments

    bwa mem \
    data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    data/fastq/mother_R1.fastq.gz \
    data/fastq/mother_R2.fastq.gz \
    > alignments/mother.sam
    ```

### 3. Alignment statistics

**Exercise:** Check out the statistics of the alignment by using `samtools flagstat`. Write the output of samtools flagstat to a file called `mother.sam.flagstat`. Do this by creating a script called `04_get_alignment_statistics.sh`, and add this script to `~/workdir/scripts/mother_only`. Find the documentation of `samtools flagstat` [here](http://www.htslib.org/doc/samtools-flagstat.html). Any duplicates in there?

??? done "Answer"

    ```sh title="04_get_alignment_statistics.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/alignments
    samtools flagstat mother.sam > mother.sam.flagstat
    ```

    This should result in:

    ```
    133477 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    317 + 0 supplementary
    0 + 0 duplicates
    132892 + 0 mapped (99.56% : N/A)
    133160 + 0 paired in sequencing
    66580 + 0 read1
    66580 + 0 read2
    131470 + 0 properly paired (98.73% : N/A)
    131990 + 0 with itself and mate mapped

    585 + 0 singletons (0.44% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
    ```

    No duplicates were found (`0 + 0 duplicates`). The aligner doesn't automatically flag duplicates. This needs to be done after the alignment.

### 4. Sorting and compression

Many downstream analyses require a coordinate sorted alignment file. Now, your alignment file is in the same order as the fastq file. You can coordinate sort an alignment file with `samtools sort`. You can find the documentation [here](http://www.htslib.org/doc/samtools-sort.html). 

**Exercise**: Sort the alignment file according to coordinate. In order to do this, create a script called `05_sort_alignment.sh` (in `~/workdir/scripts/mother_only`). 

??? done "Answer"
    ```sh title="05_sort_alignment.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results

    samtools sort -o alignments/mother.sorted.sam alignments/mother.sam 
    ```

!!! tip "Tip: `samtools sort` and `samtools view` can write to stdout"
    Like `bwa mem`, `samtools sort` and `samtools view` can write its output to stdout. This means that you need to redirect your output to a file with `>` or use the the output option `-o`.

The command `samtools view` is very versatile. It takes an alignment file and writes a filtered or processed alignment to the output. You can for example use it to compress your SAM file into a BAM file. Let's start with that.

**Exercise**: compress our SAM file into a BAM file and include the header in the output. For this, use the `-b` and `-h` options. Perform the calculation from a script called `06_compress_alignment.sh` (in `~/workdir/scripts/mother_only`).  Find the required documentation [here](http://www.htslib.org/doc/samtools-view.html). How much was the disk space reduced by compressing the file?

??? done "Answer"
    ```sh title="06_compress_alignment.sh"
    #!/usr/bin/env bash
    
    cd ~/workdir/results

    samtools view -bh alignments/mother.sorted.sam > alignments/mother.bam
    ```
    By using `ls -lh`, you can find out that `mother.sorted.sam` has a size of 55 Mb, while `mother.bam` is only 16 Mb.  

### 5. Recap

You have now performed:

- alignment
- sorting 
- compression
- flag statistics 

On the sample `mother`. Your scripts directory should look like this:

```
scripts
├── 01_download_course_data.sh
├── 02_create_bwa_index.sh
└── mother_only
    ├── 03_alignment.sh
    ├── 04_get_alignment_statistics.sh
    ├── 05_sort_alignment.sh
    └── 06_compress_alignment.sh

1 directory, 6 files
```