# 05. Data filtering and annotations
Despite using the `--filter-threshold` flag when running `modkit`, the output bedMethyl file still contains a lot of noise. Preliminary results in the bedMethyl file can be further filtered down using python. The following python script can be used to filter, organize, and export the methylation data in a `.tsv` file format (note that `pandas` is required to run this script):
* `filter-pileup.py`

Run the python script:
```sh
python filter-pileup.py pileup.bed
```

To annotate the modification sites at a single-nucleotide resolution, `snippy-vcf_to_tab` can be used. However, this requires a `.vcf` file input therefore, the `.tsv` file output from `filter-pileup.py` must be formatted and converted into a `.vcf` file format. This can be done with the python script below:

* `format-vcf.py`

Run the python script:
```sh 
python format-vcf.py 6mA.tsv
```

Run `snippy-vcf_to_tab` to annotate the `.vcf` files:
```sh
snippy-vcf_to_tab --ref {ref.fna} --gff {ref.gff} --vcf 6mA.vcf > 6mA.tab
```
Using the `.tab` file exported by `snippy-vcf_to_tab`, the annotations can be merged with the `.tsv` file from `filter-pileup.py` using the the script below:
* `annotate-tsv.py`

Run the python script:
```sh
python annotate-tsv.py 6mA.tsv snippy/6mA.tab
```