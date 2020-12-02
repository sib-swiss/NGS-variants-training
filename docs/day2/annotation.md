
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/VariantAnnotation.pdf){: .md-button }

## Exercises

To use the human genome is a reference, we have downloaded the database with:

!!! warning "No need to download, it's already downloaded for you"

```sh
# don't run this. It's already downloaded for you
snpEff download -v GRCh38.99
```

You can run snpEff like so:

```sh
snpEff -Xmx4g \
-v \
GRCh38.99 \
variants/trio.filtered.vcf > annotation/trio.filtered.snpeff.vcf
```

**Exercise:** Run the command, and check out the html file (`snpEff_summary.html`). Try to answer these questions:

**A.** How many effects were calculated?

**B.** How many variants are in the vcf?

**C.** Why is this different?

**D.** How many effects are there resulting in a missense mutation?

??? done "Answer"
    A. There were 10,357 effects calculated.

    B. There are only 556 variants in the vcf.

    C. This means that there are multiple effects per variant. snpEff calculates effects for each splice variant, and therefore the number of effects are a multitude of the number of variants.

    D. Two effects result in a missense mutation.

You can (quick and dirty) query the annotation vcf (`trio.filtered.snpeff.vcf`) for the missense mutation with `grep`.

**Exercise:** Find the variant causing the missense mutation (the line contains the string `MISSENSE`). And answer the following questions:

??? hint "Hint"
    ```sh
    grep MISSENSE annotation/trio.filtered.snpeff.vcf
    ```

**A.** How are the SNP annotations stored in the vcf?

**B.** What are the genotypes of the individuals?

**C.** Which amino acid change does it cause?

??? done "Answer"

    Find the line with the missense mutation like this:

    ```sh
    grep MISSENSE annotation/trio.filtered.snpeff.vcf
    ```

    This results in (long line, scroll to the right to see more):

    ```
    chr20   10049540        .       T       A       220.29  PASS    AC=1;AF=0.167;AN=6;BaseQRankSum=-6.040e-01;DP=85;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=8.16;ReadPosRankSum=0.226;SOR=0.951;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|cTg/cAg|L324Q|ANKEF1|protein_coding|CODING|ENST00000378392|7)      GT:AD:DP:GQ:PL   0/0:34,0:34:99:0,102,1163  0/1:17,10:27:99:229,0,492       0/0:24,0:24:72:0,72,811
    ```

    A. SNP annotations are stored in the INFO field, starting with `EFF=`

    B. The genotypes are homozygous reference for the father and son, and heterozygous for the mother. (find the order of the samples with `grep ^#CHROM`)

    C. The triplet changes from cTg to cAg, resulting in a change from L (Leucine) to Q (Glutamine).
