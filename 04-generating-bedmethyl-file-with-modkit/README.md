# 04. Generating a bedMethyl file via modkit
Official documentation for ONT's `modkit` can be found [here](https://nanoporetech.github.io/modkit/). `modkit` can be installed according to the documentation, or it can be installed in a conda environment via the following commands:
```sh
# 1. Create and activate a new conda environment
conda create -n modkit _env
conda activate modkit_env

# 2. Install modkit via the bioconda channel
conda install -c conda-forge -c bioconda ont-modkit

# 3. Verify it was installed properly
modkit --version
```

Run `pileup` to generate a bedMethyl file from the `.mapped.bam` file. A bedMethyl file is a `.bed` file containing the potential base modification sites on a single-nucleotide resolution.
```sh
modkit pileup -t 12 barcode01.mapped.bam pileup.bed -r {ref.fna} --filter-threshold 0.75
```
***
### Additional information
The bedMethyl file does not have any headers to the columns, however, their meaning can be deciphered using the `modkit` documentation linked above. It is important to note that the 4th column is the identified for what the modification was detected, however they have their own ID nomenclature for denoting the modifications. Here, looking at the methylation modifications specifically, they are as follows:

| Modification | Code | Description |
|--------------|------|-------------|
| 6mA | a | N6-methyladenine |
| 5mC | m | 5-methylcytosine |
| 4mC | 21839 | N4-methylcytosine|

More modifications and their associated ID nomenclature are detailed in the **`dorado`** documentation [here](https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/).