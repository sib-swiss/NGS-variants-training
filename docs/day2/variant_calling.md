
## Learning outcomes

**After having completed this chapter you will be able to:**

- Perform basic calculations regarding the genotype likelihood of individual variants
- Follow `gatk` best practices workflow to perform a variant analysis by:
    - Calling variants with `gatk HaplotypeCaller`
    - Combining multiple `vcf` files into a single `vcf` file
- Perform basic operations to get statistics of a `vcf` file

## Material

The [paper](https://www.biorxiv.org/content/10.1101/201178v3) on genomic variant call format (gVCF)

[GATK best practices germline short variant workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-):

<figure>
  <img src="../../assets/images/gatk_germline.png" width="800"/>
</figure>

## Exercises

### 1. Variant calling

#### Calculating PL and GQ by hand 

Here's a function in R to calculate genotype likelihoods as described in Li H.  Bioinformatics. 2011;27:2987–93 (assuming equal base error probabilities for all reads):

```R
genotype_likelihood <- function(m,g,e,ref,alt){
  (((m-g)*e+g*(1-e))^alt * ((m-g)*(1-e)+g*e)^ref)/(m^(ref+alt))
}
```

Where:

- `m` : ploidy
- `g` : number of alternative alleles
- `e` : base error probability
- `ref` : number of reference alleles counted
- `alt` : number of alternative alleles counted

**Exercise:** In a local R session, calculate the three genotype likelihoods (for g = 0, g = 1 and g = 2) for a case where we count 22 reference alleles and 4 alternative alleles (so a coverage of 26), and base error probability of 0.01. Calculate the PL values (`-10*log10(likelihood)`) for each genotype. 

!!! note "No local R installation?"
    If you don't have access to an R installation, you can also install it in the jupyter environment by running inside the terminal:

    ```sh
    conda install r-base 
    ```

    After installation, start an interactive R session by typing `R`. 

??? done "Answer"
    ```R
    # For g = 0 (i.e. 0 reference alleles)
    -10*log10(genotype_likelihood(m = 2, g= 0, e = 0.01, ref = 22, alt = 4))
    # [1] 80.96026 
    -10*log10(genotype_likelihood(m = 2, g= 1, e = 0.01, ref = 22, alt = 4))
    # [1] 78.2678
    -10*log10(genotype_likelihood(m = 2, g= 2, e = 0.01, ref = 22, alt = 4))
    # [1] 440.1746
    ``` 

**Exercise:** What is the most likely genotype? What is the genotype quality (GQ)? Do you think we should be confident about this genotype call?

??? done "Answer" 
    The most likely genotype has the lowest PL, so where g=1 (heterozygous). GL is calculated by subtracting the lowest PL from the second lowest PL, so 80.96 - 78.27 = 2.69. 
    
    This is a low genotype quality (note that we're in the phred scale), i.e. an error probability of 0.54. This makes sense, if the genotype is heterozygous we would roughly expect to count as many reference as alternative alleles, and our example quite strongly deviates from this expectation. 

#### Calling variants with GATK

The command [`gatk HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) is the core command of `gatk`. It performs the actual variant calling.

**Exercise:** Check out the [`gatk HaplotypeCaller` documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), and find out which arguments are required.

??? done "Answer"
    Required arguments are:

    * `--input`
    * `--ouput`
    * `--reference`

**Exercise:** Make a directory `~/workdir/variants` to write the output vcf. After that, run `gatk HaplotypeCaller` with required options on the recalibrated alignment file of the mother (`bqsr/mother.recal.bam`). We'll focus on a small region, so add `--intervals chr20:10018000-10220000`.

??? done "Answer"
    ```sh
    cd ~/workdir
    mkdir variants

    gatk HaplotypeCaller \
    --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    --input bqsr/mother.recal.bam \
    --output variants/mother.HC.vcf \
    --intervals chr20:10018000-10220000
    ```

**Exercise:** You can get the number of records in a vcf with piping the output of `grep -v '^#'` to `wc -l`. Get the number of variants in the vcf.

??? done "Answer"
    ```sh
    grep -v '^#' variants/mother.HC.vcf | wc -l
    ```

    Shows you that there are 411 variants in there.

You can get some more statistics with `gatk VariantsToTable`. The output can be used to easily query things in `R` or MS Excel.

Here's an example:

```sh
gatk VariantsToTable \
--variant variants/mother.HC.vcf \
--fields CHROM -F POS -F TYPE -GF GT \
--output variants/mother.HC.table
```

**Exercise:** Run the command and have a look at the first few records (use e.g. `head` or `less`). After that, report the number of SNPs and INDELs.

??? done "Answer"
    You can get the number of SNPs with:

    ```sh
    grep -c "SNP" variants/mother.HC.table
    ```

    which will give 326

    And the number of INDELs with:

    ```sh
    grep -c "INDEL" variants/mother.HC.table
    ```

    that outputs 84

    A more fancy way to this would be:

    ```sh
    cut -f 3 variants/mother.HC.table | tail -n +2 | sort | uniq -c
    ```

    Giving:

    ```
    84 INDEL
    1 MIXED
    326 SNP
    ```

We will do the variant calling on all three samples. Later we want to combine the variant calls. For efficient merging of vcfs, we will need to output the variants as a GVCF. To do that, we will use the option `--emit-ref-confidence GVCF`. Also, we'll visualise the haplotype phasing with IGV in the next section. For that we'll need a phased bam. You can get this output with the argument `--bam-output`.

**Exercise:** Run `gatk HaplotypeCaller` for mother, father and son by using a loop, and by using the arguments in the previous exercise. On top of that add the arguments `--emit-ref-confidence GVCF` and `--bamoutput <phased.bam>`.

??? done "Answer"
    ```sh
    cd ~/workdir

    for sample in mother father son
    do
      gatk HaplotypeCaller \
      --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
      --input bqsr/$sample.recal.bam \
      --output variants/$sample.HC.g.vcf \
      --bam-output variants/$sample.phased.bam \
      --intervals chr20:10018000-10220000 \
      --emit-ref-confidence GVCF
    done

    ```

### 2. Combining GVCFs

Now that we have all three GVCFs of the mother, father and son, we can combine them into a database. We do this because it enables us to later add GVCFs (with the option `--genomicsdb-update-workspace-path`), and to efficiently combine them into a single vcf.

You can generate a GenomicsDB on our three samples like this:

```sh
cd ~/workdir

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
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--variant gendb://genomicsdb \
--intervals chr20:10018000-10220000 \
--output variants/trio.vcf
```

**Exercise:** Run this command to generate the combined vcf.