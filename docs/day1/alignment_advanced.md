
## Learning outcomes

**After having completed this chapter you will be able to:**

- Use `samtools` to mark duplicates from an alignment file
- Use `samtools` to add readgroups to an alignment file
- Use a for loop in bash to perform the same operation on a range of files
- Use `samtools` in a pipe to efficiently do multiple operations on an alignment file in a single command

## Material

* `samtools` [documentation](http://www.htslib.org/doc/samtools.html)

## Exercises


### 1. Adding readgroups

During several steps of variant calling `gatk` uses read group information. For each read, this gives information on the sequencing platform, the library, the lane and of course the sample. Have a look at the description of the different levels of read group information `gatk` uses [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups). 

**Exercise:** The documentation mentions several read group fields that are used by `gatk`. Have a look at the `fastq` header. Does that give you the information that is required? Do we have that information for our sample? Can you specify it for our sample?

!!! Hint
    You can have a look at the first few entries in the `fastq` file with:

    ```sh
    zcat mother_R1.fastq.gz | head
    ```

??? done "Answer"
    Most of the information you should now based on the experimental design, the rest you can find in the [fastq header](https://en.wikipedia.org/wiki/FASTQ_format):

    - `PL`: the platform. Should be quite obvious; you usually you have this information. For us, this would be `ILLUMINA`
    - `SM`: the sample. All alignments that have reads coming from the same individual should have the same identifier in this field. For us, this would be `mother`. 
    - `LB`: library identifier. Molecular duplicates only exist within a library. If a single library was sequenced on multiple lanes, it is important to track this information. In our case, we have sequenced only one library, so you can specify it with e.g. `lib1`. 
    - `PU`: platform unit. This field is used to identify the sequencing lane. The documentation tells us we should specify it as `[FLOWCELL].[LANE].[SAMPLE BARCODE]`. The header of the first entry in our fastq file looks like this: `@H0164ALXX140820:2:1101:2136:40460/1`. Where the flowcell ID is `H0164` and the lane `2`. This formatting is specific to Broad Genomic Services pipelines, and not very common nowadays. Here the sample barcode is added to the flowcell ID, and is therefore specified as ALXX140820. We can therefore specify it with `H0164.2.ALXX140820`. 
    - `ID`: read group id. If you don't have specific information on the flowcell and lane (specified with `PU`), you can use this field to specify a unique unit that is used for e.g. base quality score recalibration. This often a combination of a flow cell identifier and a lane. In our case this could be `H0164.2`

!!! Note
    More modern output of an Illumina sequencer looks e.g. like this (example on [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format)):

    ```
    @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    ```

    Here, e.g. the `PU` field would be `FC706VJ.2.ATCACG`


**Exercise:** Have a look at the [documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-) of `AddOrReplaceReadGroups`. Specify the required arguments, and run the command. Do this from a script called `B05_add_readgroups.sh` (in `~/workdir/scripts/B-mother_only`).

??? done "Answer"
    We can use the answers of the previous exercise, and use them in the command:

    ```sh title="B05_add_readgroups.sh"
    #!/usr/bin/env bash

    cd ~/project/results

    gatk AddOrReplaceReadGroups \
    --INPUT alignments/mother.bam \
    --OUTPUT alignments/mother.rg.bam \
    --RGLB lib1 \
    --RGPU H0164.2.ALXX140820 \
    --RGPL ILLUMINA \
    --RGSM mother \
    --RGID H0164.2
    ```

**Exercise:** Compare the header and first alignments of `mother.bam` and `mother.rg.bam`. Notice any differences?

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
    @RG     ID:H0164.2      LB:lib1 PL:ILLUMINA     SM:mother       PU:H0164.2.ALXX140820
    ```

    In the alignment records, a tag was added at the very end of each line: `RG:Z:H0164.2`. Note that all fields (`LB`, `PU`, etc.) are related to `ID`. So for each read only `ID` is specified and all other fields can be deducted from that. 

### 2. Mark duplicates

Now that we have specified read groups, we can mark the duplicates with `gatk MarkDuplicates`. 

**Exercise:** Have a look at the [documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-), and run `gatk MarkDuplicates` with the three required arguments. Do this from a script called `B06_mark_duplicates.sh` (in `~/workdir/scripts/B-mother_only`). 

??? done "Answer"
    ```sh title="B06_mark_duplicates.sh"
    #!/usr/bin/env bash

    cd ~/project/results

    gatk MarkDuplicates \
    --INPUT alignments/mother.rg.bam \
    --OUTPUT alignments/mother.rg.md.bam \
    --METRICS_FILE alignments/marked_dup_metrics_mother.txt 
    ```

**Exercise:** Run `samtools flagstat` on the alignment file with marked duplicates, and write the output to a file called `mother.rg.md.bam.flagstat`. Create a script called `B07_get_alignment_stats_after_md.sh` (in `~/workdir/scripts/B-mother_only`). How many reads were marked as duplicate?

??? done "Answer"
    ```sh title="B07_get_alignment_stats_after_md.sh"
    #!/usr/bin/env bash

    cd ~/project/results/alignments

    samtools flagstat mother.rg.md.bam > mother.rg.md.bam.flagstat
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

### 3. Indexing

To look up specific alignments, it is convenient to have your alignment file indexed. An indexing can be compared to a kind of 'phonebook' of your sequence alignment file. Indexing can be done with `samtools` as well, but it first needs to be sorted on coordinate (i.e. the alignment location). You can do it like this:

```sh
samtools index <bam file>
```

**Exercise**: Create a script called `B08_index_alignment.sh` (in `~/workdir/scripts/B-mother_only`) to perform the alignment.

??? done "Answer"
    ```sh title="B08_index_alignment.sh"
    #!/usr/bin/env bash

    cd ~/project/results/alignments/

    samtools index mother.rg.md.bam
    ```


### 4. Recap - mother only

Now we have performed now the following steps on the sample `mother`:

- Read alignment
- Adding readgroups
- Marking of duplicates
- Indexing

Your `scripts` directory should look like this:

```
├── A-prepare_references
│   ├── A01_download_course_data.sh
│   └── A02_create_bwa_index.sh
├── B-mother_only
│   ├── B01_alignment.sh
│   ├── B02_get_alignment_statistics.sh
│   ├── B03_sort_alignment.sh
│   ├── B04_compress_alignment.sh
│   ├── B05_add_readgroups.sh
│   ├── B06_mark_duplicates.sh
│   ├── B07_get_alignment_stats_after_md.sh
│   └── B08_index_alignment.sh
└── C-all_samples

```

### 5. Apply it on all three samples with pipes and loops

We now apply these steps to all three samples. In order to do that, we combine the alignment with the sorting and compression in one command. We can do that with piping the output of `bwa` to `samtools sort` and `samtools view`, like this:

```sh
SAMPLE="mother"

bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
data/fastq/"$SAMPLE"_R1.fastq.gz \
data/fastq/"$SAMPLE"_R2.fastq.gz \
| samtools sort \
| samtools view -bh > results/alignments/"$SAMPLE".bam
```

**Exercise**: Make a directory in the scripts directory `C-all_samples` (so `~/workdir/scripts/C-all_samples`). In here, create a script called `C01_alignment_sorting_compression.sh`. Within that script use the above snippet to make a loop that performs the alignment, sorting and compression for all three samples (i.e. `mother`, `father` and `son`).

??? done "Answer"

    Your `scripts` directory should look like:

    ```
    scripts
    ├── A-prepare_references
    │   ├── A01_download_course_data.sh
    │   └── A02_create_bwa_index.sh
    ├── B-mother_only
    │   ├── B01_alignment.sh
    │   ├── B02_get_alignment_statistics.sh
    │   ├── B03_sort_alignment.sh
    │   ├── B04_compress_alignment.sh
    │   ├── B05_add_readgroups.sh
    │   ├── B06_mark_duplicates.sh
    │   ├── B07_get_alignment_stats_after_md.sh
    │   └── B08_index_alignment.sh
    └── C-all_samples
        └── C01_alignment_sorting_compression.sh

    ```

    And the script: 

    ```sh title="C01_alignment_sorting_compression.sh"
    #!/usr/bin/env bash

    cd ~/project

    for SAMPLE in mother father son
    do
        bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
        data/fastq/"$SAMPLE"_R1.fastq.gz \
        data/fastq/"$SAMPLE"_R2.fastq.gz \
        | samtools sort \
        | samtools view -bh > results/alignments/"$SAMPLE".bam
    done
    ```

Now we continue with adding the readgroups. For each sample, we have to add specific information to the different readgroup fields. We can do that by looping over a tab delimited file with sample-specific information in each row. Let's create that tab-delimited file. 

**Exercise** Generate a tab-delimited file called `sample_rg_fields.txt` and store it in `~/workdir/results/`. In this file, each line should represent a sample (mother, father and son), and you specify the `SM`, `LB`, `PU` and `ID` fields. E.g., the first line (for 'mother') would look like:

```
mother	lib1	H0164.2.ALXX140820	H0164.2
```

!!! warning
    Make sure to add a newline (++enter++) at the end of the file. Otherwise a loop will stop at the second-last line.

??? done "Answer"
    Your file should look like this:

    ```
    mother	lib1	H0164.2.ALXX140820	H0164.2
    father	lib2	H0164.3.ALXX140820	H0164.3
    son	lib3	H0164.6.ALXX140820	H0164.6
    ```

**Exercise** Generate a script called `C02_add_readgroups.sh` (in `~/workdir/scripts/C-all_samples`) to loop over the tab-delimited file (have a look at the last exercise in [Setup](../server_login#loops)), and add the correct readgroups to the bam file of each sample with `gatk AddOrReplaceReadGroups`. 

??? hint
    
    Try to just print the variables from a loop in order to check to see whether the loop performs according to your expectation. E.g.:

    ```sh
    cd ~/project/results

    cat sample_rg_fields.txt | while read SAMPLE LB PU ID
    do
        echo $SAMPLE $LB $PU $ID
    done
    ```

??? done "Answer"
    
    ```sh title="C02_add_readgroups.sh"
    #!/usr/bin/env bash

    cd ~/project/results

    cat sample_rg_fields.txt | while read SAMPLE LB PU ID
    do
        gatk AddOrReplaceReadGroups \
        --INPUT alignments/"$SAMPLE".bam \
        --OUTPUT alignments/"$SAMPLE".rg.bam \
        --RGLB "$LB" \
        --RGPU "$PU" \
        --RGPL ILLUMINA \
        --RGSM "$SAMPLE" \
        --RGID "$ID"
    done 
    ```

As final step, we will mark the duplicates and perform the indexing for the three samples. 

**Exercise:** Generate two scripts called `C03_mark_duplicates.sh` and `C04_index_alignments.sh`, in which you loop over the sample names and  perform the respective calculations. You can use `B06_mark_duplicates.sh` and `B08_index_alignment.sh` as a template. 

??? done "Answer"
    ```sh title="C03_mark_duplicates.sh"
    #!/usr/bin/env bash

    cd ~/project/results
    
    for SAMPLE in mother father son
    do
        gatk MarkDuplicates \
        --INPUT alignments/"$SAMPLE".rg.bam \
        --OUTPUT alignments/"$SAMPLE".rg.md.bam \
        --METRICS_FILE alignments/marked_dup_metrics_"$SAMPLE".txt 
    done

    ```

    ```sh title="C04_index_alignment.sh"
    #!/usr/bin/env bash

    cd ~/project/results

    for SAMPLE in mother father son
    do
        samtools index alignments/"$SAMPLE".rg.md.bam
    done  
    ```

### 6. Recap - all samples

Now, we have performed on all three samples:

- Alignment
- Sorting
- Compression
- Adding readgroups
- Marking of duplicates
- Indexing

Your scripts directory should look like this:

```
scripts
├── A-prepare_references
│   ├── A01_download_course_data.sh
│   └── A02_create_bwa_index.sh
├── B-mother_only
│   ├── B01_alignment.sh
│   ├── B02_get_alignment_statistics.sh
│   ├── B03_sort_alignment.sh
│   ├── B04_compress_alignment.sh
│   ├── B05_add_readgroups.sh
│   ├── B06_mark_duplicates.sh
│   ├── B07_get_alignment_stats_after_md.sh
│   └── B08_index_alignment.sh
└── C-all_samples
    ├── C01_alignment_sorting_compression.sh
    ├── C02_add_readgroups.sh
    ├── C03_mark_duplicates.sh
    └── C04_index_alignment.sh
```