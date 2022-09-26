## Learning outcomes

**After having completed this chapter you will be able to:**

* Understand the importance of reproducibility
* Apply some basic rules to support reproducibilty in computational research

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/reproducible_research.pdf){: .md-button }

### Some good practices for reproducibility

During the exercise you will be guided to adhere to the following basic principles for reproducibility:

1. **Execute the commands from a script** in order to be able to trace back your steps
2. **Number scripts** based on their order of execution (e.g. `01_download_reads.sh`)
3. Give your scripts a **descriptive and active name**, e.g. `06_build_bowtie_index.sh`
4. Make your scripts **specific**, i.e. do not combine many different commands in the same script
5. Refer to **directories and variables** at the beginning of the script

!!! tip "Keep re-evaluating your code structure"
    If you start a project it can be difficult to know what kind of analyses you are going to run, and how they interrelate. While working on a project, therefore re-evaluate readability of your project structure, and do not hesitate to change script numbering, names or contents. 

By adhering to these simple principles it will be relatively straightforward to re-do your analysis steps only based on the scripts, and will get you started to adhere to the [Ten Simple Rules for Reproducible Computational Research](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003285). 

## Exercises

Throughout the exercises today and tomorrow we will work on three different 'subprojects':

- Preparing references by indexing
- Alignment and variant calling on one sample ('mother')
- Alignment, variant calling and filtering on all samples

We store the scripts required for these subprojects in different subdirectories of `~/workdir/scripts` named:

- `A-prepare_references`
- `B-mother_only`
- `C-all_samples`

You can already create these directories now with:

```sh
cd ~/workdir/scripts/

mkdir -p \
A-prepare_references \
B-mother_only \
C-all_samples
```

By the end of day 2 `~/workdir/scripts` should look (something) like this:

```
scripts
├── A_prepare_references
│   ├── A01_download_course_data.sh
│   ├── A02_create_bwa_index.sh
│   ├── A03_create_vcf_indices.sh
│   └── A04_create_fasta_index.sh
├── B_mother_only
│   ├── B01_alignment.sh
│   ├── B02_get_alignment_statistics.sh
│   ├── B03_sort_alignment.sh
│   ├── B04_compress_alignment.sh
│   ├── B05_add_readgroups.sh
│   ├── B06_mark_duplicates.sh
│   ├── B07_get_alignment_stats_after_md.sh
│   ├── B08_index_alignment.sh
│   ├── B09_perform_bqsr.sh
│   ├── B10_run_haplotypecaller.sh
│   └── B11_variants_to_table.sh
└── C_all_samples
    ├── C01_alignment_sorting_compression.sh
    ├── C02_add_readgroups.sh
    ├── C03_mark_duplicates.sh
    ├── C04_index_alignment.sh
    ├── C05_perform_bqsr.sh
    ├── C06_run_haplotypecaller.sh
    ├── C07_create_genomicsdb.sh
    ├── C08_genotype_gvcfs.sh
    ├── C09_select_SNPs.sh
    ├── C10_select_INDELs.sh
    ├── C11_filter_SNPs.sh
    ├── C12_filter_INDELs.sh
    ├── C13_merge_filtered.sh
    ├── C14_extract_mother_only.sh
    ├── C15_evaluate_concordance.sh
    ├── C16_extract_mother_before_filtering.sh
    └── C17_evaluate_concordance_before_filtering.sh
```