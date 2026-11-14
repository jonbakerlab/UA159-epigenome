# 01. Installing Dorado and basecalling models
**Oxford Nanopore Technology's Dorado basecalling software is computationally intensive and should be run on a high-performace computing cluster.**

Download the appropriate dorado package according to ONT's [Dorado Documentation](https://software-docs.nanoporetech.com/dorado/latest/). For this project, `dorado-1.4.0-linux-x64` was used:
```sh
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz" -o dorado-1.4.0-linux-x64.tar.gz
tar -xzf dorado-1.4.0-linux-x64.tar.gz
dorado-1.4.0-linux-x64/bin/dorado --version
```
`dorado download --list` provides all models available to download for the current version of dorado being used. Use `dorado download --model` to download the super high accuracy basecalling models, for example:
```sh
path/to/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
```
Basecalling model:
* `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`

Modification models:
* `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1`
* `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v`

## Additional information
When running the basecalling for personal use, **download the most recent versions of dorado and basecalling models** for the most up-to-date and accurate basecalling.  The new versions and updates can be found in ONT's [Dorado Documentation](https://software-docs.nanoporetech.com/dorado/latest/) along with many other uses of dorado.