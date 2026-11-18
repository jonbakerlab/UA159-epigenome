# 05. Data filtering and annotations
**Placeholder variables are denoted by the `{}` and must be replaced by an appropriate file path or value.**


1. Further filter, organize, and export the bedMethyl file into a `.tsv` file format using...
  * `filter-pileup.py`

Run the python script:
```sh
python filter-pileup.py {bedMethyl.bed}
```
This code separates the bedMethyl file based on modification types. For example, if 4mC, 5mC, and 6mA modifications were basecalled, it would output three separate files:
* `4mC.tsv`
* `5mC.tsv`
* `6mA.tsv`

Steps 2-4 below, therefore, need to be done for each modification type **individually**.

2. Format and convert the `.tsv` file into a `snippy`-compatible `.vcf` format using...
   * `format-vcf.py`

Run the python script:
```sh 
python format-vcf.py {.tsv}
```

3. Generate an annotated `.tab` file by inputting the `.vcf` in `snippy-vcf_to_tab`

```sh
snippy-vcf_to_tab --ref {ref.fna} --gff {ref.gff} --vcf {.vcf} > {.tab}
```

4. Merge the `.tab` annotation to the bedMethyl-derived `.tsv` file using...
  * `annotate-tsv.py`

Run the python script:
```sh
python annotate-tsv.py {.tab} snippy/{.tab}
```

## Additional information
The script `filter-pileup.py` is coded to **only** identify 4mC, 5mC, and 6mA modifications. However by modulating lines 43 and 44, it can accommodate any modification basecalled, provided the proper SAM code (see the [04 Additional information section](04-generating-bedmethyl-file-with-modkit/) for for more information of SAM codes). For example, if 5-Hydroxymethylcytosine modifications were also basecalled, lines 43 and 44 would be adjusted as such:
```py
modifications = ["a", "m", "21839", "h"]
modification_names = {
    "a": "6mA",
    "m": "5mC",
    "21839": "4mC",
    "h": "5hmC"
}
```
