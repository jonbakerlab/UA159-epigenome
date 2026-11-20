# 06. Identifying modification motifs via modkit
**Placeholder variables are bracketed (`{}`) and must be replaced by an appropriate file path or value.**

`modkit find-motifs` is a function of modkit that uses a machine learning algorithm to tease out different motifs in a reference genome given modifications found in the bedMethyl file. This function exports a `.tsv` containing all identified motifs. **The depth of the sequencing heavily impacts the likelihood of a motif being identified and reported**. This command can be executed as follows (thread count `-t` can be adjusted as needed):

```sh
modkit find-motifs -t {threads} --in-bedmethyl {bedMethyl.bed} --ref {ref.fna} -o modkit_motifs.tsv
```

## Additional information
The interpretation of these tables, as well as how these values are determined, are described in the modkit documentation [here](https://nanoporetech.github.io/modkit/intro_find_motifs.html).

More recent versions of `modkit` have a suit of motif-related functions in `modkit motif`. Within that suit of functions is `modkt motif search`, which is functionally identical to `modkit find-motifs`. Other functions within the suit are described in the modkit documentation [here](https://nanoporetech.github.io/modkit/intro_motif.html).

### Motif confidence on a per read basis
The caveat with using `modkit find-motifs` is that it does **not** look at the probabilities of a modification (MM/ML values) on a **per read basis**, but rather what is compiled and filtered in the bedMethyl file from `modkit pileup`. If there is a motif that is less confident or wants to be further investigated, `modkit extract full` is a function that provides analytics of **each** MM/ML value in **each read** in the `.mapped.bam`, void of any compiling and filtering. This can help reduce any likelihood of false positive or false negative motifs. This analysis pipeline can be done as follows:

1. Use `seqkit` to retroactively identify the location of all instances of the motif sequence in the genome:
```sh
seqkit locate --degenerate -p {motif} {ref.fna} > {motif_sites.tsv}
```
The flag `--degenerate` allows for this function to identify IUPAC nomenclature.

2. Convert the `.tsv` file to a `.bed` file containing the location of the modified bases in the genome:
```sh
awk 'NR>1 {print $1"\t"$5+0"\t"$5+1}' {motif_sites.tsv} > {modification_sites.bed}
```
In `'NR>1 {print $1"\t"$5+0"\t"$5+1}'`, the `+0` and `+1` must be adjusted to accommodate the position of the modification site in the motif sequence.

3. Run `modkit extract full` to generate a file (`{extracts.tsv}`) containing analytics of the MM/ML values:
```sh
modkit extract full --include-bed {modification_sites.bed} --threads {threads} {.mapped.bam} {extracts.tsv}
```
The column `mod_qual` provides the value for the quality or confidence of the modification on a base.

4. Generate a histogram of the `mod_qual` values using `modqual_stats.py`:
```sh
python modqual_stats.py {extracts.tsv}
```
This will provide a visual distribution of the quality of the modifications of **each** MM/ML value in **each read**.
