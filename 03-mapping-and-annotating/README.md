# 03. Mapping and annotating the `.bam` files
Convert the `.bam` to a `.fast`. The `-T MM,ML` flag retains all modification base types and probabiliy information when converting from `.bam` to `.fastq`:
```sh
samtools fastq {.bam} -T MM,ML > {.fastq}
```
Map the `.fastq` to a `.mapped.bam` file (thread count can be adjusted): 
```sh
minimap2 -t 12 --secondary=no -ax map-ont -y {ref.fna} {.fastq} | \
samtools view -b | \
samtools sort -@ 12 - {.mapped.bam}
```
Index the `.mapped.bam` file:
```
samtools index {.mapped.bam}
```
