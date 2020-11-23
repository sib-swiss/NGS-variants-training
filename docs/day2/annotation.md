
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/VariantAnnotation.pdf){: .md-button }

## Exercises

```sh
snpEff -Xmx4g \
-v -o gatk \
GRCh38.99 \
variants/mother.trio.filtered.vcf > annotation/mother.trio.filtered.snpeff.vcf
```
