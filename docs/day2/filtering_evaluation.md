
## Learning outcomes

**After having completed this chapter you will be able to:**

- Explain why using Variant Quality Score Recalibration (VQSR) for filtering variants can outperform hard filtering 
- Perform hard filtering on both SNPs and INDELs separately by using `gatk SelectVariants` in combination with `gatk VariantFiltration`
- Perform concordance between called variants and a truth set and evaluate performance of a variant calling workflow

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/filtering_evaluation.pdf){: .md-button }

## Exercises

### 1. Hard filtering

The developers of `gatk` strongly advise to do [Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-) for filtering SNPs and INDELs. However, this is not always possible. For example, in the case of limited data availability and/or in the case you are working with non-model organisms and/or in the case you are a bit lazy and okay with a number of false positives.

Our dataset is too small to apply VQSR. We will therefore do hard filtering instead.

#### Splitting SNPs and INDELs

First, filtering thresholds are usually different for SNPs and INDELs. Therefore, we will split `trio.vcf` into two vcfs, one containg only SNPs, and one containing only INDELs. You can extract all the SNP records in our trio vcf like this:

```sh
cd ~/workdir

gatk SelectVariants \
--variant variants/trio.vcf \
--select-type-to-include SNP \
--output variants/trio.SNP.vcf
```

**Exercise:** Check out the [documentation of `gatk SelectVariants`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants), and:

* Figure out what you'll need to write at `--select-type-to-include` if you want to select only INDELS.
* Make a script (named `C09_select_SNPs.sh`) to generate a vcf with only the SNPs
* Make a script (named `C10_select_INDELs.sh`) to generate a second vcf with only the INDELs from `trio.vcf`.

??? done "Answer"
    You will need to fill in `INDEL` at `--select-type-to-include` to filter for INDELs.

    To get the SNPs you can run the command above:

    ```sh title="C09_select_SNPs.sh"
    #!/usr/bin/env bash

    cd ~/workdir

    gatk SelectVariants \
    --variant results/variants/trio.vcf \
    --select-type-to-include SNP \
    --output results/variants/trio.SNP.vcf
    ```
    
    To get the INDELs you'll need to change `--select-type-to-include` to `INDEL`:

    ```sh title="C10_select_INDELs.sh"
    #!/usr/bin/env bash

    cd ~/workdir

    gatk SelectVariants \
    --variant results/variants/trio.vcf \
    --select-type-to-include INDEL \
    --output results/variants/trio.INDEL.vcf

    ```

#### Filtering SNPs

The command `gatk VariantFiltration` enables you to filter for both the INFO field (per variant) and FORMAT field (per genotype). For now we're only interested in filtering variants. Below you can find the command to hard-filter the SNP variants on some sensible thresholds (that are explained [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)).

```sh
gatk VariantFiltration \
--variant variants/trio.SNP.vcf \
--filter-expression "QD < 2.0"              --filter-name "QD2" \
--filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
--filter-expression "SOR > 3.0"             --filter-name "SOR3" \
--filter-expression "FS > 60.0"             --filter-name "FS60" \
--filter-expression "MQ < 40.0"             --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output variants/trio.SNP.filtered.vcf
```

**Exercise:** Run the filtering command above in a script called `C11_filter_SNPs.sh`. Did it affect the number of records in the vcf?

!!! hint
    You can check out the number of records in a vcf with:

    ```sh
    grep -v "^#" <variants.vcf> | wc -l
    ```

??? done "Answer"
    Your script:

    ```sh title="C11_filter_SNPs.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/variants

    gatk VariantFiltration \
    --variant trio.SNP.vcf \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --output trio.SNP.filtered.vcf
    ```

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
--filter-expression "QD < 2.0"                  --filter-name "QD2" \
--filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
--filter-expression "FS > 200.0"                --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
--output variants/trio.INDEL.filtered.vcf
```

**Exercise:** Run the command from a script called `C12_filter_INDELs.sh` and figure out how many variants are filtered out.

!!! hint
    You can use this command from the answer to the previous exercise:

    ```sh
    grep -v "^#" <variants.vcf> | cut -f 7 | sort | uniq -c
    ```

    to see how many INDELs were filtered out.

??? done "Answer"

    Your script:

    ```sh title="C12_filter_INDELs.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/variants

    gatk VariantFiltration \
    --variant trio.INDEL.vcf \
    --filter-expression "QD < 2.0"                  --filter-name "QD2" \
    --filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
    --filter-expression "FS > 200.0"                --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
    --output trio.INDEL.filtered.vcf
    ```

    And check out the contents of the `FILTER` column: 

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
--INPUT <input1.vcf> \
--INPUT <input2.vcf> \
--OUTPUT <merged.vcf>
```

**Exercise:** Run this command from a script called `C13_merge_filtered.sh` to merge the vcfs (`trio.SNP.filtered.vcf` and `trio.INDEL.filtered.vcf`).

??? done "Answer"

    ```sh title="C13_merged_filtered.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/variants

    gatk MergeVcfs \
    --INPUT trio.SNP.filtered.vcf \
    --INPUT trio.INDEL.filtered.vcf \
    --OUTPUT trio.filtered.vcf

    ```


### 2. Evaluation by concordance

For this region we have a highly curated truth set for the mother available. It originates from the [Illumina Platinum truth set](https://www.illumina.com/platinumgenomes.html). You can find it at `data/variants/NA12878.vcf.gz`

To check how well we did, we'd first need to extract a vcf with only the information of the mother.

**Exercise:** Generate a script called `C14_extract_mother_only.sh` to extract variants that have at least one alternative allele in the mother from `variants/trio.filtered.vcf`. Use `gatk SelectVariants` with the arguments:

* `--sample-name mother`
* `--exclude-non-variants`
* `--remove-unused-alternates`

In addition to the required arguments.

??? done "Answer"
    ```sh title="C14_extract_mother_only.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/variants

    gatk SelectVariants \
    --variant trio.filtered.vcf \
    --sample-name mother \
    --exclude-non-variants \
    --remove-unused-alternates \
    --output mother.trio.filtered.vcf
    ```

**Exercise:**

A. How many variants are in the extracted vcf? How many of those are filtered out?

B. Compare our vcf with the curated truth set with the command below from a script called `C15_evaluate_concordance.sh`. How many SNPs didn't we detect?

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

    Your script to evaluate the concordance:

    ```sh title="C15_evaluate_concordance.sh"
    #!/usr/bin/env bash

    cd ~/workdir

    gatk Concordance \
    --evaluation results/variants/mother.trio.filtered.vcf \
    --truth data/variants/NA12878.vcf.gz \
    --intervals chr20:10018000-10220000 \
    --summary results/variants/concordance.mother.trio.filtered

    ```

    Check out the output with `cat`:

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

**Exercise:** Check out the concordance of the mother with the truth set before filtering. Do this by generating two scripts:

- `C16_extract_mother_before_filtering.sh`: to run `gatk SelectVariants` in order to get only variants from the mother from the unfiltered `trio.vcf`. 
- `C17_evaluate_concordance_before_filtering.sh`: to run `gatk Concordance` on the selected variants. 

 Did filtering improve the recall or precision?

!!! note
    We did the filtering on `trio.vcf`, therefore, you first have to extract the records that only apply to the mother by using `gatk SelectVariants`.

    Also note that `trio.vcf` contains records other than SNPs and INDELs. Use `--select-type-to-include` to select only SNPs and INDELs.

??? done "Answer"

    First select only SNPs and INDELs from the mother from the unfiltered vcf:

    ```sh title="C16_extract_mother_before_filtering.sh"
    #!/usr/bin/env bash

    cd ~/workdir/results/variants

    gatk SelectVariants \
    --variant trio.vcf \
    --sample-name mother \
    --exclude-non-variants \
    --remove-unused-alternates \
    --select-type-to-include INDEL \
    --select-type-to-include SNP \
    --output mother.trio.vcf

    ```

    Get the concordance with the truth set:

    ```sh title="C17_evaluate_concordance_before_filtering.sh"
    #!/usr/bin/env bash

    cd ~/workdir

    gatk Concordance \
    --evaluation results/variants/mother.trio.vcf \
    --truth data/variants/NA12878.vcf.gz \
    --intervals chr20:10018000-10220000 \
    --summary results/variants/concordance.mother.trio


    ```

    Which gives:

    ```
    type    TP      FP      FN      RECALL  PRECISION
    SNP     319     7       9       0.973   0.979
    INDEL   63      20      6       0.913   0.759
    ```

    The precision for SNPs is slightly lower. Due to filtering, we removed two false positives.
