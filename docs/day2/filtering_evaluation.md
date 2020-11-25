
## Exercises

### 1. Hard filtering

The developers of `gatk` strongly advise to do [Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-) for filtering SNPs and INDELs. However, this is not always possible. For example, in the case of limited data availability and/or in the case you are working with non-model organisms and/or in the case you are a bit lazy and okay with a number of false positives.

Our dataset is too small to apply VQSR. We will therefore do hard filtering instead.

#### Splitting SNPs and INDELs

First, filtering thresholds are usually different for SNPs and INDELs. You can extract all the SNP records in our trio vcf like this:

```sh
cd ~/workdir

gatk SelectVariants \
--variant variants/trio.vcf \
--select-type SNP \
--output variants/trio.SNP.vcf
```

**Exercise:** Check out the [documentation of `gatk SelectVariants`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants), and:

* Figure out what you'll need to fill in at `--select-type` if you want to select only INDELS.
* Generate a vcf with only the SNPs and a second vcf with only the INDELs from `trio.vcf`.

??? done "Answer"
    You will need to fill in `INDEL` at `--select-type` to filter for INDELs.

    To get the SNPs you can run the command above. To get the INDELs you'll need to change `--select-type` to `INDEL`:

    ```sh
    gatk SelectVariants \
    --variant variants/trio.vcf \
    --select-type INDEL \
    --output variants/trio.INDEL.vcf
    ```

#### Filtering SNPs

The command `gatk VariantFiltration` enables you to filter for both the INFO field (per variant) and FORMAT field (per genotype). For now we're only interested in filtering variants. Below you can find the command to hard-filter the SNP variants on some sensible thresholds (that are explained [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)).

```sh
gatk VariantFiltration \
--variant variants/trio.SNP.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output variants/trio.SNP.filtered.vcf
```

**Exercise:** Run the filtering command above. Did it affect the number of records in the vcf?

!!! hint
    You can check out the number of records in a vcf with:

    ```sh
    grep -v "^#" <variants.vcf> | wc -l
    ```

??? done "Answer"
    There are no differences in the number of records:

    ```sh
    grep -v "^#" variants/trio.SNP.vcf | wc -l
    ```

    and

    ```sh
    grep -v "^#" variants/trio.SNP.filtered.vcf | wc -l
    ```

    both give 446.

    However, there are SNPs filtered out, by changing the `FILTER` column. You can check the number of records with PASS by:

    ```sh
    grep -v "^#" variants/trio.SNP.filtered.vcf | cut -f 7 | sort | uniq -c
    ```

    Giving:

    ```
    441 PASS
    2 QD2;SOR3
    3 SOR3
    ```

#### Filtering INDELs

A command with [sensible parameters](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering) to do a first iteration of hard filtering the INDELs would be:

```sh
gatk VariantFiltration \
--variant variants/trio.INDEL.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--output variants/trio.INDEL.filtered.vcf
```

**Exercise:** Run the command and figure out how many variants are filtered out.

!!! hint
    You can use this command from the answer to the previous exercise:

    ```sh
    grep -v "^#" <variants.vcf> | cut -f 7 | sort | uniq -c
    ```

    to see how many INDELs were filtered out.

??? done "Answer"
    ```sh
    grep -v "^#" variants/trio.INDEL.filtered.vcf | cut -f 7 | sort | uniq -c
    ```

    gives:

    ```
    110 PASS
    ```

    So no variants are filtered out.

#### Merging filtered SNPs and INDELs

Now that we have filtered the INDELs and SNPs separately, we can merge them again with this command:

```sh
gatk MergeVcfs \
--INPUT variants/trio.SNP.filtered.vcf \
--INPUT variants/trio.INDEL.filtered.vcf \
--OUTPUT variants/trio.filtered.vcf
```

**Exercise:** Run the command to merge the vcfs.

### 2. Evaluation by concordance

For this region we have a highly curated truth set for the mother available. It originates from the [Illumina Platinum truth set](https://www.illumina.com/platinumgenomes.html). You can find it at `variants/NA12878.vcf.gz`

To check how well we did, we'd first need to extract a vcf with only the information of the mother.

**Exercise:** To extract variants that have at least one alternative allele in the mother from `variants/trio.filtered.vcf`, use `gatk SelectVariants` with the arguments:

* `--sample-name mother`
* `--exclude-non-variants`
* `--remove-unused-alternates`

In addition to the required arguments.

??? done "Answer"
    ```sh
    gatk SelectVariants \
    --variant variants/trio.filtered.vcf \
    --sample-name mother \
    --exclude-non-variants \
    --remove-unused-alternates \
    --output variants/mother.trio.filtered.vcf
    ```

**Exercise:**

A. How many variants are in `mother.trio.filtered.vcf`? How many of those are filtered out?

B. Compare our vcf with the curated truth set with the command below. How many SNPs didn't we detect?

```sh
gatk Concordance \
--evaluation variants/mother.trio.filtered.vcf \
--truth data/variants/NA12878.vcf.gz \
--intervals chr20:10018000-10220000 \
--summary variants/concordance.mother.trio.filtered
```

??? done "Answer"

    To get the number of records per FILTER, we run:

    ```sh
    grep -v "^#" variants/mother.trio.filtered.vcf | cut -f 7 | sort | uniq -c
    ```

    gives:

    ```
    407 PASS
    2 SOR3
    ```

    So two records were filtered out, based on the Symmetric Odds Ratio (issues with strand bias).

    Check out the output of `gatk Concordance` with `cat`:

    ```sh
    cat variants/concordance.mother.trio.filtered
    ```

    gives:

    ```
    type    TP      FP      FN      RECALL  PRECISION
    SNP     319     5       9       0.973   0.985
    INDEL   63      20      6       0.913   0.759
    ```

    Showing that there were 9 false negatives, i.e. SNPs we didn't detect.

!!! tip "Recall & precision"
    More info on the definition of recall and precision on this [wikipedia page](https://en.wikipedia.org/wiki/Precision_and_recall)

**Exercise:** Check out the concordance of the mother with the truth set before filtering. Did filtering improve the recall or precision?

!!! note
    We did the filtering on `trio.vcf`, therefore, you first have to extract the records that only apply to the mother by using `gatk SelectVariants`.

    Also note that `trio.vcf` contains records other than SNPs and INDELs. Use `--select-type` to select only SNPs and INDELs.

??? done "Answer"

    First select only SNPs and INDELs from the mother from the unfiltered vcf:

    ```sh
    gatk SelectVariants \
    --variant variants/trio.vcf \
    --sample-name mother \
    --exclude-non-variants \
    --remove-unused-alternates \
    --select-type INDEL \
    --select-type SNP \
    --output variants/mother.trio.vcf
    ```

    Get the concordance with the truth set:

    ```sh
    gatk Concordance \
    --evaluation variants/mother.trio.vcf \
    --truth data/variants/NA12878.vcf.gz \
    --intervals chr20:10018000-10220000 \
    --summary variants/concordance.mother.trio
    ```

    Which gives:

    ```
    type    TP      FP      FN      RECALL  PRECISION
    SNP     319     7       9       0.973   0.979
    INDEL   63      20      6       0.913   0.759
    ```

    The precision for SNPs is slightly lower. Due to filtering, we removed two false positives.
