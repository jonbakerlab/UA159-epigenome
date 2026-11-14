# 07. Constructing a consensus logo
Consensus logos visualize the relative abundance of bases surrounding all of the modification sites of interest. The python script below takes the position of each modification position within a `.tsv` file filtered from the original bedMethyl file and, using a provided `.fasta` file, will pull the nucleotide at the given modification position as well as the bases +/- 12 positions upstream and downstream of the modification. This is iterated for each position present in the `.tsv`, and will generate a consensus logo compiling all sequence windows surrounding the modified bases. `biopython` python package is required. This can downloaded by the following command:

```sh
pip install biopython
```
The python script to generate a consensus logo is as follows:
* `consensus_logo.py`

Run the python script:
```sh
python consensus_logo.py annotated_table_6mA.tsv {ref.fna} {strain_name}
```