
## Learning outcomes

**After having completed this chapter you will be able to:**

- Use `samtools` to mark duplicates from an alignment file
- Use `samtools` to add readgroups to an alignment file
- Use a for loop in bash to perform the same operation on a range of files
- Use `samtools` in a pipe to efficiently do multiple operations on an alignment file in a single command

## Material

* `samtools` [documentation](http://www.htslib.org/doc/samtools.html)

## Exercises

### 1. Marking duplicates

For variant analysis, it's important to mark reads that possibly originated from PCR duplication. We can do that with `samtools markdup`. However, we can not directly run that on our `.sam` file nor on our compressed `.bam` file.

**Exercise:** Which samtools commands would we need to run to mark duplicates?

!!! hint
    More info on this at the [`samtools markdup` documentation](http://www.htslib.org/doc/samtools-markdup.html#EXAMPLES).

    Also: the reads are already collated (i.e. forward and reverse are grouped) after the alignment, so no need to run `samtools collate`.

??? done "Answer"
    The commands would be:

    * `samtools fixmate` (with the `-m` option)
    * `samtools sort`
    * `samtools markdup`

**Exercise:** Run the three commands that are required to mark duplicates.

??? done "Answer"
    ```sh
    cd ~/workdir/alignment

    samtools fixmate -m mother.bam mother.fixmate.bam
    samtools sort -o mother.positionsort.bam mother.fixmate.bam
    samtools markdup mother.positionsort.bam mother.markdup.bam
    ```

** Exercise:** Run `samtools flagstat` on the alignment file with marked duplicates. How many reads were marked as duplicate?

??? done "Answer"
    ```sh
    samtools flagstat mother.markdup.bam
    ```

    Gives:

    ```
    133477 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    317 + 0 supplementary
    17329 + 0 duplicates
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
    Which tells us that 17329 reads were marked as duplicate.

### 2. Adding read groups

For variant analysis, it's important to know which read came from which sample. Right now, that's easy. All reads come from one individual. But this can become less trivial if you are combining samples. Therefore we add a tag to each read specifying its origin.

You can add a readgroup to your marked alignment file like this:

```sh
samtools addreplacerg \
-r ID:mother \
-r SM:mother \
-r PL:ILLUMINA \
-o mother.markdup.rg.bam \
mother.markdup.bam
```

This command modifies the sam header and read tags.

**Exercise:** Run the `samtools addreplacerg` command and compare the header and first alignments of `mother.markdup.bam` and `mother.markdup.rg.bam`. Notice any differences?

!!! hint
    You can view the header with

    ```
    samtools view -H <alignment.bam>
    ```

    And the first few alignments with

    ```
    samtools view <alignment.bam> | head
    ```

??? done "Answer"
    Compared to the header of `mother.markdup.bam`, the header of `mother.markdup.rg.bam` contains an extra line starting with `@RG`:

    ```
    @RG     ID:mother       SM:mother       PL:ILLUMINA
    ```

    In the alignment records, a tag was added at the very end of each line: `RG:Z:mother`.


### 3. Indexing

To look up specific alignments, it is convenient to have your alignment file indexed. An indexing can be compared to a kind of 'phonebook' of your sequence alignment file. Indexing can be done with `samtools` as well, but it first needs to be sorted on coordinate (i.e. the alignment location). You can do it like this:

```sh
samtools index mother.markdup.rg.bam
```

### 4. Piping and looping

Samtools can quite easily be used in a UNIX pipeline. This has the advantage that you don't need to write many intermediate files. However, the developers have not been very consistent with managing input and output (I'm sure they had their reasons). To use samtools in a pipe, the input argument needs to replaced with a `-`. Also, some commands do not write by default to stdout, but to a specified file (this is the case for e.g. `samtools fixmate` and `samtools markdup`). In that case, also the output argument should be replaced with `-`.

!!! example
    For the command `samtools addreplacerg` the samtools documentation provides the following synopsis:

    ```
    samtools addreplacerg [-r rg-line | -R rg-ID] [-m mode] [-l level] [-o out.bam] in.bam
    ```

    Meaning that it requires `in.bam`, and can write to `out.bam` if option `-o` is provided. By default, it writes to stdout. So, if you pipe to `samtools addreplacerg` you would only need to replace the input file with a `-`:

    ```sh
    some_command | samtools addreplacerg [options] - > output.sam
    ```

    For the command `samtools fixmate`, the samtools documentation provides this synopsis:

    ```
    samtools fixmate [-rpcm] [-O format] in.nameSrt.bam out.bam
    ```

    Meaning that it requires both the input file `in.nameSrt.bam` and the output file `out.bam`. So, if you pipe to `samtools fixmate` and you want to write to stdout (so piping from), you'll need to replace both the input and output with a `-`:

    ```sh
    some_command | samtools fixmate [options] - - > output.sam
    ```

    The most frequently used samtools commands don't require an input nor an output file, and therefore behave like many UNIX commands. An example of this is `samtools sort`. The synopsis is:

    ```
    samtools sort [-l level] [-m maxMem] [-o out.bam] [-O format] [-n] [-t tag] [-T tmpprefix] [-@ threads] [in.sam|in.bam|in.cram]
    ```

    Telling us that both the input file and output (with option `-o`) file are optional. If the input file is absent, it reads from stdin. So, you could use it without a `-` replacing input or output files:

    ```sh
    some_command | samtools sort > output.sam
    ```

Let's put everything we've done so far in a pipe and loop over our three samples.

The command below loops over the strings `father`, `mother` and `son`, and performs these tasks:

1. Create a variable to work on data of mother, father and son separately
2. Perform the alignment
3. Fill in the mate coordinates and sort on coordinate
4. Mark duplicates
5. Add readgroups
6. Compress the output
7. Create an index

```sh
#!/usr/bin/env bash

cd ~/workdir

for sample in mother father son
do
  bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
  data/fastq/"$sample"_R1.fastq.gz \
  data/fastq/"$sample"_R2.fastq.gz \
  | samtools fixmate -m - - \
  | samtools sort \
  | samtools markdup -s - - \
  | samtools addreplacerg -r ID:$sample -r SM:$sample -r PL:ILLUMINA - \
  | samtools view -bh > alignment/$sample.bam

  samtools index alignment/$sample.bam
done
```

**Exercise:** For each task (1-7), figure out which part of the script performs that task. After that, run it to get the alignments of all three samples.

??? done "Answer"
    Creating variables (1):

    ```sh
    for sample in mother father son
    do
      ...
    done
    ```

    Perform the alignment (2):

    ```sh
    bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    data/fastq/"$sample"_R1.fastq.gz \
    data/fastq/"$sample"_R2.fastq.gz
    ```

    Fill in the mate coordinates and sort on coordinate (3):

    ```sh
    samtools fixmate -m <INPUT> <OUTPUT> \
    | samtools sort
    ```

    Mark duplicates (4):

    ```sh
    samtools markdup -s <INPUT> <OUTPUT>
    ```

    Add readgroups (5):

    ```sh
    samtools addreplacerg -r ID:$sample -r SM:$sample -r PL:ILLUMINA <INPUT>
    ```

    Compress the output (6):

    ```sh
    samtools view -bh > alignment/$sample.bam
    ```

    Create an index (7):

    ```sh
    samtools index alignment/$sample.bam
    ```
