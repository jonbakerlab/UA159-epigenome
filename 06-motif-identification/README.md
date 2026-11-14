# 06. Identifying modification motifs via modkit
`modkit find-motifs` is a function of modkit that uses a machine learning algorithm to tease out different motifs in a reference genome given modifications found in the bedMethyl file. This function exports a `.tsv` containing all identified motifs. **The depth of the sequencing heavily impacts the likelihood of a motif being identified and reported**. This command can be executed as follows (thread count `-t` can be adjusted as needed):

```sh
modkit find-motifs -t 12 --in-bedmethyl pileup.bed --ref {ref.fna} -o modkit_motifs.tsv
```
***
### Additional information
The interpretation of these tables, as well as how these values are decided, are descrived in the modkit documentation [here](https://nanoporetech.github.io/modkit/intro_find_motifs.html).

More recent versions of `modkit` have a suit of motif-related functions in `modkit motif`. Within that suit of functions is `modkt motif search`, which is functionally identical to `modkit find-motifs`. Other functions within the suit are described in the modkit documentation [here](https://nanoporetech.github.io/modkit/intro_motif.html).