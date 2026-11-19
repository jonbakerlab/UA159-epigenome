# 10. RNAseq analysis using deseq2
**Placeholder variables are bracketed (`{}`) and must be replaced by an appropriate file path or value.**

RNAseq is a powerful tool that can investigate differentially expressed genes between different samples. This is a step-by-step process on how to take raw `.fastq` files from a sequencing run, into a compatible gene-expression-matrix format for deseq2.

## Mapping the sequencing data
**A reference `.fna` is required to map the raw `.fastq` files.**
1. Index the `.fna` file using `bwa-mem2`
```sh
bwa-mem2 index {ref.fna}
```
2. Map the `.fast` files via `minimap2`
    * Run `sbatch_sam.sh`
3. Convert the `.sam` files to `.bam` files via `samtools`
    * Run `sbatch_bam.sh`
4. Sort the `.bam` files via `samtools`
    * Run `sbatch_sort.sh`
5. Index the `.mapped.bam` files via `samtools`

**All** scripts must be edited and adjusted to accommodate appropriate usage needs and file paths. All placeholder variables are denoted by brackets `{}` (unless preceded by a `$`, i.e., `${}` as that bash language for calling a variable). Importantly, `#SBATCH --array=` must contain the range of sample barcodes ID numbers. For example: `--array=8-19`.

## Creating the feature table
**A reference `.gtf` is required to use `subreds featureCounts`.**

With all of the `.mapped.bam` files for each sample, `subreads featureCOunts` can be used to generate the gene expression matrix (featureCounts table) so be inputted and analyzed by deseq2. The function can be run as follows:

```sh
#!/bin/bash
conda activate subread_env
srun --cpus-per-task 4 featureCounts \
    -a {ref.gtf} \
    -t gene \
    -g gene_id \
    --extraAttributes gene,product,gene_biotype \
    -L \
    -O \
    -M \
    --primary \
    -T {threads} \
    -o featureCounts.txt \
    *.mapped.bam
```
With the outputted `featureCounts.txt` file, it can be formatted properly manually to (via excel or some other method of editing a `.txt` file) such that it is formatted accordingly:

| Geneid | sample1_1 | sample1_2 | sample1_3 | sample2_1 | sample2_2 | sample2_3 | sampleN | ... |
|--------|-----------|-----------|-----------|-----------|-----------|-----------|---------|-----|
| SMU_RS00005 | 18 | 14 | 20 | 83 | 35 | 40 | n | ... | 

**Note the flag `-L`, which is important as it signifies that these are long read sequences. This applies for reads that are sequenced on a nanopore, whereas, if it was an illumina sequence, the sequences are short reads.**

## Run deseq2 on R
With the `featureCounts.txt` file properly formatted, it can then be run in `R` to perform a deseq2 analysis via `deseq2-pairwise-comparison.R` To properly run the script the following lines must be edited to apply to work for the `featureCounts.txt` input file:
* **Line 15**: Input an appropriate working directory in `setwd("~/path/to/working/dir/")`
    - Ex:  `setwd("~/RNAseq/Smutans/methylation-mutants)`
* **Line 22**: Make sure the proper file name is present, otherwise the feature table will not be read.
    - Ex: `countdata <- read.table("featureCount.txt", header=TRUE, row.names=1)`
* **Line 30**: Provide the path to the reference `.gtf` file in `{ref.gtf}`
    - Ex:  `../NCBI-files/GCF_000007465.2_ASM746v2_genomic.gtf`
* **Line 75**: Input the samples/conditions and the number of replicates in the format provided.
    - Ex: 
    ```R
    (condition <- factor(c(
        rep("UA159", 3), 
        rep("dpnII", 3),
        rep("hsdM", 3)
    )))
    ```
* **Line 94**: Input the appropriate reference sample/condition in `dds$condition <- relevel(dds$condition, "{reference_condition}")`.
    - Ex: `dds$condition <- relevel(dds$condition, "UA159")`

Assuming everything is in order, this will output a multitude of files containing deseq2 results for each piecewise comparison, volcano plots, heatmaps, and much more.
