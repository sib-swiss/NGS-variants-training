
## Learning outcomes

**After having completed this chapter you will be able to:**

- Describe the aims of variant annotation
- Explain how variants are ranked in order of importance
- Explain how splice variation affects variant annotation
- Perform a variant annotation with `snpEff`
- Interpret the report generated by `snpEff`
- Explain how variant annotation can be added to a `vcf` file

## Material

Presentation will be sent to you by e-mail.

## Exercises

To use the human genome as a reference, we have downloaded the database with:

!!! warning "No need to download, it's already downloaded for you"

```sh
# don't run this. It's already downloaded for you
snpEff download -v GRCh38.99
```

You can run snpEff like so:

```sh
mkdir annotation

snpEff -Xmx4g \
-v \
-o gatk \
GRCh38.99 \
variants/trio.filtered.vcf > annotation/trio.filtered.snpeff.vcf
```

!!! note "Output `-o gatk` is deprecated for `gatk4`"
    Here, we use output `-o gatk` for readability reasons (only one effect per variant is reported). With `gatk3` you could use `gatk VariantAnnotator` with input from `snpEff`. In `gatk4` that is not supported anymore.

**Exercise:** Run the command, and check out the html file (`snpEff_summary.html`). Try to answer these questions:

**A.** How many effects were calculated?

**B.** How many variants are in the vcf?

**C.** Why is this different?

**D.** How many effects result in a missense mutation?

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

!!! note "Only one effect per SNP in the vcf"
    In the vcf we have created you can only find one effect per SNP. If you would run `snpEff` without `-o gatk`, you would get all effects per variant.

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