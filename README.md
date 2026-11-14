# UA159 epigenome

This repository accompanies our upcoming publication on the epigenome of *Streptococcus mutans* UA159. It contains workflow documentation, scripts, and configuration files for identifying and analyzing DNA methylation using Oxford Nanopore sequencing data.

## Workflow sections

The numbered directories describe the analysis from basecalling through visualization. Each section includes its own instructions and any accompanying scripts or configuration files.

1. [Installing Dorado and basecalling models](01-installing-dorado-and-basecalling-models/) - Install Dorado and the basecalling and modified-base models.
2. [Basecalling](02-basecalling/) - Basecall POD5 files with modified-base models and demultiplex samples.
3. [Mapping and annotating](03-mapping-and-annotating/) - Map reads to a reference genome while retaining modification tags.
4. [Generating a bedMethyl file with modkit](04-generating-bedmethyl-file-with-modkit/) - Summarize modification calls from mapped reads.
5. [Filtering and annotations](05-filtering-and-annotations/) - Filter modification calls and annotate their genomic locations.
6. [Motif identification](06-motif-identification/) - Identify modification-associated sequence motifs with modkit.
7. [Constructing sequence logos](07-constructing-sequence-logos/) - Visualize sequence context around modification sites.
8. [Hemi- versus fully methylated sites](08-hemi-versus-fully-methylated/) - Compare methylation at paired sites on opposite strands.
9. [Generating Circos plots](09-generating-circos-plots/) - Visualize the genomic distribution of modification sites.

Start with the individual section READMEs for commands, dependencies, and input requirements. Basecalling is computationally intensive and uses GPUs; examples and paths may need to be adapted to your data and computing environment.
