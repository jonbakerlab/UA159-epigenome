# 03. Mapping modBAM files
**Placeholder variables are denoted by the `{}` and must be replaced by the path to the appropriate file.**
1. Convert the `.bam` to a `.fast` using `samtools`: 
```sh
samtools fastq {.bam} -T MM,ML > {.fastq}
```

The `-T MM,ML` flag retains all modification base types and probability information **on a per read basis** when converting from `.bam` to `.fastq`. More information on the flags and their meaning can be found in the [Modified Basecalling](https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/#sam-tags-mm-ml) page in the Dorado documentation.

2. Map the `.fastq` to a `.mapped.bam` file using a pipeline of `minimap2` and `samtools`: 
```sh
minimap2 -t {threads} --secondary=no -ax map-ont -y {ref.fna} {.fastq} | \
samtools view -b | \
samtools sort -@ {threads} -o {.mapped.bam}
```
The `-t` and `-@` flags for `minimap2` and `samtools`, respectively, are for thread counts; these can be adjusted as needed.

3. Index the `.mapped.bam` file:
```
samtools index {.mapped.bam}
```
