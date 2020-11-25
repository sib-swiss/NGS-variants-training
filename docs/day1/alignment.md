## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_alignment.pdf){: .md-button }

## Exercises

### 1. Prepare the reference genome

Download and unpack the data files in your working directory (`~/workdir`).

```sh
cd ~/workdir
wget https://ngs-variants-training.s3.eu-central-1.amazonaws.com/ngs-variants-training.tar.gz
tar -xvf ngs-variants-training.tar.gz
```

This will create the directory `data`. Check out what's in there.

We'll use `bwa mem` for the alignment. Like all alignment software, it requires an index of the reference genome. You can make an index like this:

```sh
bwa index <reference.fa>
```

Make an index of the reference sequence of chromosome 20 of the human genome. You can find the fasta file in `~/workdir/data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa`.

??? done "Answer"
    ```sh
    cd ~/workdir/data/reference/
    bwa index Homo_sapiens.GRCh38.dna.chromosome.20.fa
    ```

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
    bwa mem \
    ~/workdir/data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    <forward_reads.fq> \
    <reverse_reads.fq> \
    > <alignment.sam>
    ```

Perform an alignment with `bwa mem` of the reads from the mother (`mother_R1.fastq` and `mother_R2.fastq`) against chromosome 20. Write the alignment file to a directory in `~/workdir` called `alignment`.

!!! note "Index prefix is the same a reference filename"
    With default values, the index of a reference for `bwa mem` is the same as the reference itself. In this case, this would be `Homo_sapiens.GRCh38.dna.chromosome.20.fa`.

??? done "Answer"
    We'll first make the alignment directory:

    ```sh
    cd ~/workdir/
    mkdir alignment
    ```
    Then, we run the alignment:

    ```sh
    bwa mem \
    data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    data/fastq/mother_R1.fastq \
    data/fastq/mother_R2.fastq \
    > alignment/mother.sam
    ```

### 2. Alignment statistics

**Exercise:** Check out the statistics of the alignment by using `samtools flagstat`. Find the documentation [here](http://www.htslib.org/doc/samtools-flagstat.html). Any duplicates in there?

??? done "Answer"
    ```sh
    cd ~/workdir/alignment
    samtools flagstat mother.sam
    ```

    Should give:

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




### 3. Compression

The command `samtools view` is very versatile. It takes an alignment file and writes a filtered or processed alignment to the output. You can for example use it to compress your SAM file into a BAM file. Let's start with that.

**Exercise**: compress our SAM file into a BAM file and include the header in the output. For this, use the `-b` and `-h` options. Find the required documentation [here](http://www.htslib.org/doc/samtools-view.html). How much was the disk space reduced by compressing the file?

!!! tip "Tip: `samtools view` writes to stdout"
    Like `bwa mem`, `samtools view` writes its output to stdout. This means that you need to redirect your output to a file with `>` or use the the output option `-o`.

??? done "Answer"
    ```sh
    samtools view -bh mother.sam > mother.bam
    ```
    By using `ls -lh`, you can find out that `mother.sam` has a size of 54 Mb, while `mother.bam` is only 20 Mb.  
