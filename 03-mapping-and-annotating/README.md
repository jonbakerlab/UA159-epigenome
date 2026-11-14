# 03. Mapping and annotating the `.bam` files
Convert the `.bam` to a `.fast`. The `-T MM,ML` flag retains all modification base types and probabiliy information when converting from `.bam` to `.fastq`:
```sh
samtools fastq barcode01.bam -T MM,ML > barcode01.fastq
```
Map the `.fastq` to a `.mapped.bam` file (thread count can be adjusted): 
```sh
minimap2 -t 12 --secondary=no -ax map-ont -y {ref.fna} barcode01.fastq | \
samtools view -b | \
samtools sort -@ 12 - barcode01.mapped.bam
```
Index the `.mapped.bam` file:
```
samtools index barcode01.mapped.bam
```