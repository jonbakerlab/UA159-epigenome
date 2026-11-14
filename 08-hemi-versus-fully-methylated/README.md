# 08. Elucidating fully- and hemi-methylated sites
Because of the complementing nature of the forward and reverse strands in DNA, a motif sequence on the forward strand is also found on the reverse strand. This leads to the possibility of a methylation site of a given motif on the forward strand having a pairing methylation site on the reverse strand motif sequence. This can be referred to as fully-methylation and hemi-methylation:
* **fully-methylation**: when there is a methylation site on both the forwards and reverse strand for a given motif.
* **hemi-methylation**: when there is a methylation site on only the forward **or** reverse strand for a given motif.

The `fully-hemi-methylation.py` script takes a +/- 12 base pair window surrounding the methylation site and filters the sites based on the motif provided (IUPAC nomenclature okay). Then it reports the total number of methylation sites, total number of methylation sites found under the motif provided, and the number of fully- and hemi-methylated sites based on if there is a methylation site present the correct number of bases away on the reverse strand.
* `fully-hemi-methylation.py`

Run the python script:
```py
python fully-hemi-methylation.py annotated_table_6mA.tsv {ref.fna} CGANNNNNNNTCY 2
```

### Additional information
This script also has the capacity to export `.tsv` files of all of the methylation sites with the extracted sequence (`df`), all of the methylation sites filtered by the motif provided (`df_motif`), or which of those remaining methylation sites are fully-methylated (`df_fully`) or hemi-methylated (`df_hemi`). To do this, simply remove the commenting symbol, `#`, next to the `.to_csv` line of code to export the `.tsv` of interest.