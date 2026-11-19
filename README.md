# UA159 epigenome

This repository accompanies our upcoming publication on the epigenome of *Streptococcus mutans* UA159. It contains workflow documentation, scripts, and configuration files for identifying and analyzing DNA methylation using Oxford Nanopore sequencing data.

## Workflow sections

The numbered directories describe the analysis from basecalling through visualization. Each section includes its own instructions and any accompanying scripts or configuration files.

1. [Installing Dorado and basecalling models](01-installing-dorado-and-basecalling-models/) - Install Dorado and the basecalling and modified-base models.
2. [Basecalling](02-basecalling/) - Basecall POD5 files with modified-base models and demultiplex samples.
3. [Mapping modBAM files](03-mapping-and-annotating/) - Map reads to a reference genome while retaining modification tags.
4. [Generating a bedMethyl file with modkit](04-generating-bedmethyl-file-with-modkit/) - Summarize modification calls from mapped reads.
5. [Filtering and annotations](05-filtering-and-annotations/) - Filter modification calls and annotate their genomic locations.
6. [Motif identification](06-motif-identification/) - Identify modification-associated sequence motifs with modkit.
7. [Constructing sequence logos](07-constructing-sequence-logos/) - Visualize sequence context around modification sites.
8. [Fully- vs hem- methylated sites](08-hemi-versus-fully-methylated/) - Compare methylation at paired sites on opposite strands.
9. [Generating Circos plots](09-generating-circos-plots/) - Visualize the genomic distribution of modification sites.
10. [RNAseq via deseq2](10-RNAseq-via-deseq2/) - Analyze differentially expressed genes between samples using deseq2.

Start with the individual section READMEs for commands, dependencies, and input requirements. Basecalling is computationally intensive and uses GPUs; examples and paths may need to be adapted to your data and computing environment.

## Software
| **Tool** | **Version** |
|----------|------------------|
| dorado | 1.4.0-linux-x64 |
| conda | 24.11.3 |
| python | 3.12.3 |
| R | 4.5.2 |
| circos | 0.69-8 |

| **Conda Packages** | **Version** |
|--------------------|------------------|
| samtools | 1.21 |
| minimap2 | 2.28-r120 |
| bwa-mem2 | 2.2.1 |
| subread featurCounts | 2.0.6 |
| modkit |  0.6.4 |
| snippy | 4.6.0 |

| **Python packages** | **Version** |
|---------------------|-------------|
| pandas | 2.2.3 |
| biopython | 1.85 |
| logomaker | 0.8.7 |
| matplotlib | 3.10.0 |

## Acknowledgements

