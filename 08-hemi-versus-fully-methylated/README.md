# 08. Elucidating fully- and hemi-methylated sites
**Placeholder variables are bracketed (`{}`) and must be replaced by an appropriate file path or value.**

Because of the complementing nature of the forward and reverse strands in DNA, a motif sequence on the forward strand is also found on the reverse strand. This leads to the possibility of a methylation site of a given motif on the forward strand having a pairing methylation site on the reverse strand motif sequence. This can be referred to as fully-methylation and hemi-methylation:
* **fully-methylation**: when there is a methylation site on both the forwards and reverse strand for a given motif.
* **hemi-methylation**: when there is a methylation site on only the forward **or** reverse strand for a given motif.

The `fully-hemi-methylation.py` script takes a +/- 12 base pair window surrounding the methylation site and filters the sites based on the motif provided (**IUPAC nomenclature okay**). Then it reports the total number of methylation sites, total number of methylation sites found under the motif provided, and the number of fully- and hemi-methylated sites based on if there is a methylation site present the correct number of bases away on the reverse strand.

A few notes on the input variables for this script:
* `{motif}` is an input variable that requires a string of text for the motif sequence of interest.
  - Ex: `GATC`
* `{mod_offset}` is an input variable that requires an **integer value** for the position of the modification site within the motif string in a **0-based indexing system**. This is also the same value as the `offset` value in the `modkit_motifs.tsv` file exported by `modkit find-motifs`.
  - Ex: For the input motif sequence `GATC`, the `{mod_offset}` for the methylated "A" is `1`.


Run the python script:
```py
python fully-hemi-methylation.py {.tsv} {ref.fna} {motif} {mod_offset}
```

### Additional information
This script also has the capacity to export additional `.tsv` files. To do this, simply remove the commenting symbol, `#`, from...
* **Line 260**: Export `methylation_sites_sequences.tsv`
  - A `.tsv` containing all of the methylation sites with the surrounding base pair sequence.
* **Line 263**: Export `motif_methylation_sites.tsv`
  - A `.tsv` containing all of the methylation sites that were associated with the provided motif sequence.
* **Line 266**: Exports `fully-methylated_sites.tsv`
  - A `.tsv` containing of the methylation sites associated with the motif which are considered fully-methylated.
* **Line 267**: Exports `hemi-methylated_sites.tsv`
  - A `.tsv` containing the methylation sites associated with the motif which are considered hemi-methylated.
