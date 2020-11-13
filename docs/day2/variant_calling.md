
## Exercises

### 1. Indexing, indexing, indexing

Many algorithms work faster, or only work with an index of their (large) input files. In that sense, `gatk` is no different from other tools. The index for a reference needs to be created in two steps:

```sh
cd ~/data/reference
samtools faidx <reference.fa>
gatk CreateSequenceDictionary --REFERENCE <reference.fa>
```

Also input vcf files need to be indexed. This will create a `.idx` file associated with the `.vcf`. You can do this like this:

```sh
gatk IndexFeatureFile --input <variants.vcf>
```

**Exercise:** Create the required `gatk` indexes for:

* The reference genome `reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa`
* A part of the dbsnp database: `variants/GCF.38.filtered.renamed.vcf`
* A part of the 1000 genomes indel golden standard: `variants/1000g_gold_standard.indels.filtered.vcf`

??? done "Answer"
    Creating the index for the reference genome:

    ```sh
    cd ~/data
    samtools faidx reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
    gatk CreateSequenceDictionary --REFERENCE reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
    ```

    Creating the indices for the vcfs:

    ```sh
    gatk IndexFeatureFile --input variants/1000g_gold_standard.indels.filtered.vcf
    gatk IndexFeatureFile --input variants/GCF.38.filtered.renamed.vcf
    ```

!!! note "Chromosome names"
    Unlike IGV, `gatk` requires equal chromosome names for all its input files and indexes, e.g. in `.fasta`, `.bam` and `.vcf` files. In general, for the human genome there are three types of chromosome names:

    * Just a number, e.g. `20`
    * Prefixed by `chr`. e.g. `chr20`
    * Refseq name, e.g. `NC_000020.11`

    Before you start the alignment, it's wise to check out what chromosome naming your input files are using, because changing chromosome names in a `.fasta` file is easier than in a `.bam` file.

    If your fasta titles are e.g. starting with a number you can add `chr` to it with `sed`:

    ```
    sed s/^>/>chr/g <reference.fasta>
    ```

    You can change chromsome names in a vcf with [`bcftools annotate`](http://samtools.github.io/bcftools/bcftools.html#annotate):

    ```sh
    bcftools annotate --rename-chrs <tab-delimited-renaming> <input.vcf>
    ```

### 2. Base Quality Score Recalibration (BQSR)

BQSR evaluates the base qualities on systematic error. It can ignore sites with known variants. BQSR helps to identify faulty base calls, and therefore reduces the chance on discovering false positive variant positions.

BQSR is done in two steps:

1. Recalibration with [`gatk BaseRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593511-BaseRecalibrator)
2. By using the output of `gatk BaseRecalibrator`, the application to the bam file with [`gatk ApplyBQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)

**Exercise:** Check out the documentation of the tools. Which options are required?

??? done "Answer"
    For `gatk BaseRecalibrator`:

    * `--reference`
    * `--input`
    * `--known-sites`
    * `--output`

    For `gatk ApplyBQSR`:

    * `--bqsr-recal-file`
    * `--input`
    * `--output`

**Exercise:** Run the two commands with the required options on `mother.bam`, with `--known-sites` `variants/1000g_gold_standard.indels.filtered.vcf` and `variants/GCF.38.filtered.renamed.vcf`.

!!! hint "Multiple inputs for same argument"
    In some cases you need to add multiple inputs (e.g. multiple `vcf` files) into the same argument (e.g. `--known-sites`). To provide multiple inputs for the same argument in `gatk`, you can use the same argument multiple times, e.g.:

    ```sh
    gatk BaseRecalibrator \
    --reference <reference.fa> \
    --input <alignment.bam> \
    --known-sites <variants1.vcf> \
    --known-sites <variants2.vcf> \
    --output <output.table>
    ```

??? done "Answer"
    ```sh
    cd ~/data

    mkdir bqsr

    gatk BaseRecalibrator \
    --reference reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    --input alignment/mother.bam \
    --known-sites variants/GCF.38.filtered.renamed.vcf \
    --known-sites variants/1000g_gold_standard.indels.filtered.vcf \
    --output bqsr/mother.recal.table

    gatk ApplyBQSR \
    --input alignment/mother.bam \
    --bqsr-recal-file bqsr/mother.recal.table \
    --output bqsr/mother.recal.bam
    ```

**Exercise:** Place these commands in a loop, that performs the BQSR for mother, father and son.

??? done "Answer"
    ```sh
    cd ~/data

    for sample in mother father son
    do
      gatk BaseRecalibrator \
      --reference reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
      --input alignment/$sample.bam \
      --known-sites variants/GCF.38.filtered.renamed.vcf \
      --known-sites variants/1000g_gold_standard.indels.filtered.vcf \
      --output bqsr/$sample.recal.table

      gatk ApplyBQSR \
      --input alignment/$sample.bam \
      --bqsr-recal-file bqsr/$sample.recal.table \
      --output bqsr/$sample.recal.bam
    done
    ```

### 3. Variant calling

The command [`gatk HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) is the core command of `gatk`. It performs the actual variant calling.

**Exercise:** Check out the [`gatk HaplotypeCaller` documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), and find out which arguments are required.

??? done "Answer"
    Required arguments are:

    * `--input`
    * `--ouput`
    * `--reference`

**Exercise:** Run `gatk HaplotypeCaller` with required options on the recalibrated alignment file of the mother (`bqsr/mother.recal.bam`). We'll focus on a small region, so add `--intervals chr20:10018000-10220000`.

??? done "Answer"
    ```sh
    gatk HaplotypeCaller \
    --reference reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    --input bqsr/mother.recal.bam \
    --output variants/mother.HC.vcf \
    --intervals chr20:10018000-10220000 \
    ```

We will do the variant calling on all three samples. Later we want to combine the variant calls. For efficient merging of vcfs, we will need to output the variants as a GVCF. To do that, we will use the option `--emit-ref-confidence GVCF`. Also, we'll visualise the haplotype phasing with IGV in the next section. For that we'll need a phased bam. You can get this output with the argument `--bam-output`.

**Exercise:** Run `gatk HaplotypeCaller` for mother, father and son by using a loop, and by using the arguments in the previous exercise. On top of that add the arguments `--emit-ref-confidence GVCF` and `--bamoutput <phased.bam>`.

??? done "Answer"
    ```sh
    for sample in mother father son
    do
      gatk HaplotypeCaller \
      --reference reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
      --input bqsr/$sample.recal.bam \
      --output variants/$sample.HC.g.vcf \
      --bam-output variants/$sample.phased.bam \
      --intervals chr20:10018000-10220000 \
      --emit-ref-confidence GVCF
    done

    ```

### 4. Visualisation

Download the following files to your local computer:

* `variants/mother.phased.bam`
* `variants/mother.HC.g.vcf`
* `bqsr/mother.recal.bam`

Launch IGV and select the human genome version hg38 as a reference.

<figure>
  <img src="../../assets/images/select_hg38.png" width="300"/>
</figure>

Load the downloaded files as tracks in igv with **File > Load From File...**, and navigate to region `chr20:10,026,397-10,026,638`.

**Exercise:** Zoom out for a bit. Not all reads are in the track of `mother.phased.bam`. What kind of reads are in there?

??? done "Answer"
    The reads supporting called variants.

Now, we'll investigate the haplotype phasing.  Go back to `chr20:10,026,397-10,026,638`.

!!! tip
    If your screen isn't huge, you can remove the track `mother.recal.bam`. Do that by right-click on the track, and click on **Remove Track**.


In the track with `mother.phased.bam`, right click on the reads and select **Group alignments by > read group**. This splits your track in two parts, one with artificial reads describing haplotypes that were taken in consideration (ArtificalHaplotypeRG), and one with original reads that support the variants.

**Exercise**: How many haplotypes were taken into consideration? How many haplotypes can you expect at maximum within a single individual?

!!! hint
    This might come as a shock, but humans are diploid.

??? done "Answer"
    Three haplotypes, as there are three artificial reads:

    <figure>
      <img src="../../assets/images/artificial_reads.png" width="500"/>
    </figure>

    Diploids can carry two haplotypes. So at least one of the three is wrong.

Now colour the reads by phase. Do that with by right clicking on the track and choose **Colour alignments by > tag**, and type in "HC" (that's the tag where the phasing is stored).

**Exercise:** Which artificial read doesn't get support from the original sequence reads? Are the alternative alleles of the two SNPs on the same haplotype (i.e. in phase)?

??? done "Answer"

    The track should look like this (colours can be different):

    <figure>
      <img src="../../assets/images/phased_bam.png" width="500"/>
    </figure>

    The reads only support the brown and blue haplotype, and not the pink one.

    The alternative alleles are coloured in IGV. For the first SNP this is the C (in blue) and for the second the T (in red). They are always in different reads, so they are in repulsion (not in phase).

### 5. Combining GVCFs

Now that we have all three GVCFs of the mother, father and son, we can combine them into a database. We do this because it enables us to later add GVCFs (with the option `--genomicsdb-update-workspace-path`), and to efficiently combine them into a single vcf.

You can generate a GenomicsDB on our three samples like this:

```sh
gatk GenomicsDBImport \
--variant variants/mother.HC.g.vcf \
--variant variants/father.HC.g.vcf \
--variant variants/son.HC.g.vcf \
--intervals chr20:10018000-10220000 \
--genomicsdb-workspace-path genomicsdb

```

**Exercise:** Run this command to generate the database.

You can retrieve the combined vcf from the database with `gatk GenotypeGVCFs`.

```sh
gatk GenotypeGVCFs \
--reference reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--variant gendb://genomicsdb \
--intervals chr20:10018000-10220000 \
--output variants/trio.vcf
```

**Exercise:** Run this command to generate the combined vcf.
