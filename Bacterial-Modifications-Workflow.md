# Bacterial Modification Workflow
This document provides a step-by-step, user friendly way of investigating DNA (or RNA) modifications from Oxford Nanopore Technology, Inc (ONT) sequencing technology. It will describe all the installation protocols for first time usage, as well as supplemental codes to allow for a very modular experience. This workflow uses a combination of ONT's `dorado` and `modkit` to identify and analyze modifications as well as some in-house analysis methods for interpreting and visualizing the results of these programs. So, without further ado, lets get started.

## Basecalling with methylation models and demultiplexing using ONT Dorado
**Due to the size of the raw POD5 files exported by the nanopore sequencing run, it is recommended to upload these to a computing cluster to facilitate the processing of these large file sizes. Make sure that all the POD5 files to be basecalled are all in the same directory.**

The ONT basecalling software `dorado` can be installed via the Oxford Nanopore Technology `dorado` [github page](https://github.com/nanoporetech/dorado). It is important to download the installer corresponding to the platform's operating system architecture. For example, the computing cluster used here operates on a `linux_x86_64` architecture, therefore `dorado-1.0.2-linux-x64` installer was downloaded. Once the archive has been downloaded, the archive can be extracted using `gunzip`. Note that downloading and installing `dorado` through this path does not innately put `dorado` in an executable path, therefore, make sure the archive is located in a logical and accessible directory such that it can be executed manually.

`dorado` runs by providing it with a downloaded basecalling model. The models available for a given version of dorado can be listed by running the `path/to/dorado download --list` command. Since this specific pipeline deals with sequenced DNA, the most recent and highest accuracy (denoted by the `sup` in the version name) should be downloaded. At the time of writing this (and with `dorado` version 1.0.2), model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` was downloaded and used. The model can simply be downloaded by executing the following line of code:

```sh
path/to/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

In addition to the basecalling model, supplemental modification models must be downloaded as well to basecall in a modification-aware way. Again, all the modification models are also listed using the `download --list` command. Its best practice to download the cognate modification models to the basecalling model, i.e., download the same accuracy (`sup`) and versions. So, for this example, models `dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1` and `dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1` were downloaded to identify the 6mA and the 4mC/5mC DNA methylation sites, respectively. Similar to above, these models can be downloaded using the `download --model` command:

```sh
path/to/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1
path/to/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1
```

With the proper basecalling models and modification models downloaded, the `dorado basecalling` program can be used with the `--modified-bases-models` flag to basecall in a modification-aware way. The positional argument for supplying the POD5 files to `dorado basecalling` requires the path to the **directory** containing the POD5 files, not the files themselves. Note that running the basecalling program with a super-high accuracy model (denoted by `sup` in the model name), in conjunction with the modification models, is a memory intensive process and therefore GPUs should be the primary source of memory supplemented with CPUs; more on this later. 

Since all of the barcodes are multiplexed in the POD5 files, a single BAM file with the modification information in it (referred to as modBAM) is produced. The function `dorado demux` can be used to demultiplex the modBAM file, producing a modBAM file for each barcode used/identified in a given run. The name of the kit used for the library prep and sequencing must be provided to ensure the barcodes IDs are identified properly; the name is supplied to the `--kit-name` flag, in this case it would be `NQK-NBD114-24`. The code below describes both the `dorado basecalling` and `dorado demux` process, as well as supplemental code rename and organize files.

`sbatch-modBAM-demux.sh`
```sh
#!/bin/bash

#SBATCH --job-name=dorado
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=3
#SBATCH --mem-per-gpu=20G
#SBATCH --cpus-per-gpu=4
#SBATCH --time=14:00:00
#SBATCH --export=ALL

absolute_path="/absolute/path/"
pod5_path="to/working/dir/pod5"
working_dir="to/working/dir"

mkdir ${absolute_path}/${working_dir}/dorado-output/

srun ${absolute_path}/dorado-1.0.2-linux-x64/bin/dorado basecaller \
     ${absolute_path}/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
     ${absolute_path}/${pod5_path} \
     --modified-bases-models ${absolute_path}/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1,${absolute_path}/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 \
     --modified-bases-batchsize 128 \
     --min-qscore 9 \
     --batchsize 128 \
     > ${absolute_path}/${working_dir}/dorado-output/modBAM.bam

srun ${absolute_path}/shared_resources/dorado-1.0.2-linux-x64/bin/dorado demux \
   ${absolute_path}/${working_dir}/dorado-output/modBAM.bam \
   --kit-name SQK-NBD114-24 \
   --output-dir ${absolute_path}/${working_dir}/dorado-output/demux

for f in "${absolute_path}/${working_dir}/dorado-output/demux/"*_barcode*.bam "${absolute_path}/${working_dir}/dorado-output/demux/"*_unclassified.bam; do
  newname=$(echo "$f" | sed -E 's/.*_(barcode[0-9]+|unclassified)\.bam/\1.bam/')
  mv "$f" "${absolute_path}/${working_dir}/dorado-output/demux/${newname}"
done
```
Following up on the conversation of using GPUs rather than CPUs, `dorado` can struggle with memory allocation when running the program. Meaning that if too much memory is provided (i.e. to many Gb of GPU power), it can over use the memory quickly and run into an error. Through trial and error, the parameters above seem to a good standard for minimizing total run-time while not running into an error. However, if the total size of the POD5 files supplied is very large (e.g. human genome sequences), the general rule of thumb would be to either decrease the total memory available or lowering the batch size via the `--batch-size` flag. In addition, if the run does run into an error, the flag `--resume-from` can be used with the argument being the path to the incomplete modBAM file. Note that if this is used, the output file **must** have a different name than the argument provided for `--resume-from`, e.g. `modBAM-2.bam`. This can be implimented as shown below:

```sh
srun ${absolute_path}/dorado-1.0.2-linux-x64/bin/dorado basecaller \
     ${absolute_path}/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
     ${absolute_path}/${pod5_path} \
     --modified-bases-models ${absolute_path}/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1,${absolute_path}/shared_resources/dorado-1.0.2-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 \
     --modified-bases-batchsize 128 \
     --min-qscore 9 \
     --batchsize 128 \
     --resume-from ${absolute_path}/${working_dir}/dorado-output/modBAM.bam \
     > ${absolute_path}/${working_dir}/dorado-output/modBAM-2.bam
```
## Generating a bedMethyl file using ONT Modkit
This workflow is inspired and derived from the rki-mf1/ont-methylation pipeline. The codes below are verbatum the codes used under the hood in rki-mf1/ont-methylation, however using these codes separated from the original pipeline allows for a more modular experience by selecting filtering threshold parameters, as well as being applied to other modification analyses, for example RNA modifications, rather than being exclusive to DNA methylation.

Before running `modkit` for the first time, it needs to be installed. There are two ways to do this: one uses `conda` to install it (less involved) while the other uses `rust` to install it (more involved). When initially running this workflow, `modkit` was installed using `rust`; it wasn't until later when it was discovered that it can be installed using `conda`. Therefore, when running this workflow for the first time, attempt to install `modkit` into a conda environment and see if it works. If it does, great, if not, the supplemental codes for installing `modkit` via `rust` will be provided below. To download `modkit` using `conda`, use the following steps:

```sh
conda create -n modkit_env
conda activate modkit_env
conda install -c bioconda meme
conda install -c nanoporetech modkit
```

If the above method does not work when trying to run `modkit`, the code below shows how to download and install modkit using `rust`:
```sh
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Build modkit using cargo
git clone https://github.com/nanoporetech/modkit.git
cd modkit
cargo install --path

# Install supplemental python packages
pip install pandas
pip install biopython
```

Because `modkit` requires a mapped modBAM file (and supplemental BAI file), the first order of business is to map the modBAM file(s) exported by the `dorado`, as, in the code presented above, no mapping occured during the basecalling process. To do this, the modBAM must be converted to a FASTQ which can be done using `samtools fastq`. It is critical, however, to use the flag `-T MM,ML` when running this, because this retains all of the methylation information and probabilities within the modBAM. From the FASTQ, the standard workflow of producing a mapped BAM file can be used. These steps are shown below:

```sh
samtools fastq modBAM.bam -T MM,ML > modBam.fastq
minimap2 -t 10 --secondary=no -ax map-ont -y .fna modBAM.fastq | \
samtools view -b | \
samtools sort -@ 10 -o modBAM.mapped.bam

samtools index -@ 10 modBAM.mapped.bam
```
The parameters set for the argument `-ax` in `minimap2` are dependent on what's being sequenced (gDNA, cDNA, direct RNA) so make sure to use `minimap2 --help` to determine what's best for the given data. In this case gDNA was sequenced and therefore `map-ont` will suffice.

With the mapped modBAM file, `modkit` can be run. Once again this is a fairly computationally intensive program and is best used on a computing cluster. After `modkit` is properly installed, it can be run with the following line of code:

```sh
modkit pileup -t 14 modBAM.mapped.bam pileup.bed -r .fna --filter-threshold 0.75
```
Here, the `--filter-threshold` value was selected based on the value used in the original rki-mf1/ont-methylation pipeline. The justification for that value is described in the GitHub page. `modkit` produces a BED file containing the potential methylation sites known as a bedMethyl file. If viewed, the columns in the bedMethyl file have no headers, however the meaning of each column and their values can be deciphered using the `modkit` documentation. It is important to note that the 4th column is the identifier for what modification was detected, however they have their own ID nomenclature for denoting the modifications. A few examples include:

* 5-Methylcytosine = `m`
* 6-Methyladenine = `a`
* N(4)-Methylcytosine = `21839` (based on the CHEBI index no.)

Finally, the size of the bedMethyl file can range from 5 GB to 100 GB depending on the size of the genome (i.e.,  prokaryote versus eukaryote genomes). Therefore, the remainder of the data analysis could be done locally (and therefore interactively via Jupyter notebooks) or continued on the computing cluster.

## Data filtering and annotations
Despite using the `--filter-threshold` flag when running `modkit`, the output bedMethyl file still contains a lot of noise. So, pulling similar methodologies from rki-mf1/ont-methylation, the preliminary results in the bedMethyl file can be further filtered down using python. The following python script can be used to filter, organize, and export the methylation data in a TSV format (note that `pandas` is required to run this script):

`filter-pileup.py`
```py
import os
import pandas as pd
import argparse

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Filter the bedMethyl file exported by modkit and create individual TSV files for each modification")
parser.add_argument("input1", help="Provide the path to the bedMethyl files, e.g. pileup.bed")
args = parser.parse_args()

# %% Import and format the bedMethyl file into a DataFrame
bedMethyl_df = pd.read_csv(args.input1, sep='\t', header=None)

# Name the columns based on the modkit documentation
bedMethyl_df.columns = [
    "Contig", "Position", "drop0", "Modification", "drop1", "Strand", "drop2", "drop3", "drop4", "Valid_coverage",  
    "Fraction_modified", "Modified_bases", "Unmodified_bases", "Other_mod_bases", "drop5", "Modification_below_threshold", "Other_bases", "drop6"
]

# Drop all the columns containing 'drop'
bedMethyl_df = bedMethyl_df[[c for c in bedMethyl_df.columns if not c.startswith("drop")]]

# Calculate a new parameters "Total_coverage" and "Percent_modified" (a la rki-mf1/ont-methylation protocol which will be used downstream as a filtering parameter
bedMethyl_df["Total_coverage"] = (
    bedMethyl_df["Modified_bases"]
    + bedMethyl_df["Unmodified_bases"]
    + bedMethyl_df["Other_mod_bases"]
    + bedMethyl_df["Modification_below_threshold"]
    + bedMethyl_df["Other_bases"]
)

bedMethyl_df["Percent_modified"] = bedMethyl_df["Modified_bases"] / bedMethyl_df["Total_coverage"]

# Adjust the values in the "Position" column from a 0-based index to a 1-based index
bedMethyl_df['Position'] = bedMethyl_df['Position'] + 1

# %% Filter rows from the DataFrame based on the "Total_coverage" and "Percent_modified" values and export the results
total_coverage = 10
percent_modified = 0.50

# Create a dictionary for the methylation modifications and their assiciated IDs
modifications = ["a", "m", "21839"]
modification_names = {
    "a": "6mA",
    "m": "5mC",
    "21839": "4mC"
}

# Filter the bedMethyl DataFrame as well as update the modification IDs to their actual DNA modification
filtered_table = bedMethyl_df[(bedMethyl_df.Total_coverage >= total_coverage) & (bedMethyl_df.Percent_modified >= percent_modified)].reset_index(drop=True)
for modification in modifications:
    for row_index in range(0, filtered_table.shape[0]):
        if filtered_table.loc[row_index, "Modification"] == modification:
            filtered_table.loc[row_index, "Modification"] = modification_names[modification]

# Create a new directory to hold all the filtered bedMethyl files and export them
os.makedirs("modification_tables", exist_ok=True)
filtered_table.to_csv(os.path.join("modification_tables", "filtered_table.tsv"), sep='\t', index=False)
for modification in modifications:
    modification_table = filtered_table[(filtered_table.Modification == modification_names[modification])]
    modification_table.to_csv(os.path.join("modification_tables", modification_names[modification] + ".tsv"), sep='\t', index=False)
```
A few notes about this script:
* The position values in the "Position" column of the bedMethyl DataFrame is on a 0-based indexing system. +1 was added to each of these to convert it to a 1-based indexing system due to conventional genome annotations not having a "0th" position.
* The filtering parameters can be easily modulated and changed based on what is desired; simply change the values of the "total_coverage" and "percent_modified" variables.
* The benefit of doing this independent of the rki-mf1/ont-methylation pipeline is that this can be customized and changed to accomodate a wider range of modifications that may be identified in RNA, for example. To do this, simply adjust the "modifications" list and "modfication_names" dictionary with the new modifications and their associated ID.

**The remainder of this section is completely independent of the rki-mf1/ONT-methylation pipeline**
The exported files from this filtering process contain all the modification sites deemed significant and their associated position on a single-nucleotide resultion. The next step is to annotate and determine what gene these modifications reside in and this can be done using the `vcf_to_tab` function from `snippy`. Snippy can easly be installed via the following:

```sh
conda create -n snippy_env python=3.10
conda activate snippy_env
conda install -c conda-forge -c bioconda -c defaults snippy
```

In order to use `snippy-vcf_to_tab`, the input file (in this case the filtered modkit tables) needs to be in a VCF format. This can easily be done by making placeholder columns in the TSV files and changing the file extension to `.vcf`. The python script below does just that:

`format-vcf.py`
```py
import os
import pandas as pd
import argparse

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Convert the filtered modkit table TSV files (derived from the bedMethyl file) into a VCF compatible format.")
parser.add_argument("input1", help="Provide the path to a filtered modkit table TSV file.")
args = parser.parse_args()

# %% Import the TSV file into a DataFrames
vcf_format = pd.read_csv(args.input1, sep='\t', usecols=['Contig', 'Position'])

# rename the first two columns according to the format of a VCF
vcf_format.columns.values[0] = '#CHROM'
vcf_format.columns.values[1] = 'POS'

# Add empty columns as placeholders to format the VCF file
for column in ['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']:
    vcf_format[column] = None

# Make a new directory (if not made already) and export the VCF DataFrame into it
os.makedirs('snippy', exist_ok=True)
vcf_format.to_csv(os.path.join('snippy', os.path.splitext(os.path.basename(args.input1))[0] + '.vcf'), sep='\t', index=False)
```

Now, with the formatted VCF file, it can be run through `snippy-vcf_to_tab`. Note that this program specifically requires a FASTA file (FA or FNA files also acceptable) as well as a GFF file to run.

```sh
snippy-vcf_to_tab \
 --ref .fna \
 --gff .gff \
 --vcf .vcf > .tab
```

The exported TAB file now contains the annotations for each modification position. Now, the goal is to amend those annotations to the original modification table TSV file. This can be done using the following python script:


`annotate-tsv.py`
```py
import os
import pandas as pd
import argparse

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Amend the annotations in the TAB file produced by 'snippy-vcf_to_tab' to the corresponding filtered bedMethyl TSV file")
parser.add_argument("input1", help="Path to the filtered bedMethyl TSV file of interest")
parser.add_argument("input2", help="Path to the corresponding TAB file generated by 'snippy-vcf_to_tab'")
args = parser.parse_args()

# %% Import the two input arguments into DataFrames
tsv_df = pd.read_csv(args.input1, sep='\t')
tab_df = pd.read_csv(args.input2, sep='\t')

# Add additional columns to the TSV to amend the annotations
for column in ['Locus_tag', 'Gene', 'Product']:
    tsv_df[column] = None

# %% Amend the annotations from the TAB file to the TSV file via a series of nested for-loops and if-loops

# for-loop to cycle through the rows in the TSV file DataFrame
for i in range(0, tsv_df.shape[0]):

    # Position of the i-th row
    temp_position = tsv_df.loc[i, 'Position']

    # for loop to cycle through the rows of the annotated TAB DataFrame
    for j in range(0, tab_df.shape[0]):

        # If the position is found in the j-th row of the TAB DataFrame...
        if temp_position == tab_df.loc[j, 'POS']:

            # ...then pull the information from the TAB DataFrame and tabulate it in the TSV DataFrame
            tsv_df.loc[i, 'Locus_tag'] = tab_df.loc[j, 'LOCUS_TAG']
            tsv_df.loc[i, 'Gene'] = tab_df.loc[j, 'GENE']
            tsv_df.loc[i, 'Product'] = tab_df.loc[j, 'PRODUCT']

# %% Export the annotated TSV file to the current directory
tsv_df.to_csv('annotated_table_' + os.path.splitext(os.path.basename(args.input1))[0] + '.tsv', sep='\t', index=False)
```

**Note that this current iteration of this code uses a sequence of nested for loops and if loops. So, while this is manageable with smaller file sizes, larger file sizes may take a while due to the computational intensity of this method. A better method would be to use pandas integrated function `pd.match` which does the same thing as the nested loops, however its much more memory efficient and does it within one line of code. This will be implemented in the future.**

This gives a good starting point for downstream analyses. The section below contains supplemental code for various methods to interpret these results as well as graphical visualizations.

## Supplemental analysis codes
### Venn diagram
A potential method of comparison is generating a Venn diagram to visually look at the overlapping modification sites between two strains (e.g., *S. mutans* UA159 vs *dpnII*) and as well as creating a list of the sites exclusively in one strain, exculively in the other, and the overlapping sites. This can be done easily in python with `pandas` in conjunction with `matplotlib`, which can be installed using the following:

```py 
pip install matplotlib matplotlib-venn
```

With this package installed, a Venn diagram can be created using the following script:

`venn.py`
```py
import os
import pandas as pd
import argparse
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Generate a Venn diagram of the modification sites between two samples as well as export tables of each group, i.e., in group A, in group B, and in groups A and B.")
parser.add_argument("input1", help="Path to the Sample 1 annotated TSV file.")
parser.add_argument("input2", help="Strain name for Sample 1, e.g., UA159")
parser.add_argument("input3", help="Path to the Sample 2 annotated TSV file.")
parser.add_argument("input4", help="Strain name for Sample 2, e.g., dpnII")
args = parser.parse_args()

# %% Make a Venn diagram to show the overlapping modification sites

# Import the two TSV files into DataFrames
tsv1_df = pd.read_csv(args.input1, sep='\t')
tsv2_df = pd.read_csv(args.input3, sep='\t')

# Convert 'Position columns to sets
set1 = set(tsv1_df['Position']) # Set from sample 1
set2 = set(tsv2_df['Position']) # Set from sample 2

# Use venn2 to create and export the Venn diagram
venn2([set1, set2], set_labels=(args.input2, args.input4))
os.makedirs('venn_diagram', exist_ok=True)
plt.savefig(os.path.join('venn_diagram', "venn_diagram.pdf"), format='pdf', bbox_inches='tight')

# %% Using the position sets above, generate a new DataFrame for each subgroup

# Compute subsets
only1_set = set1 - set2 # modification positions only in sample 1
only2_set = set2 - set1 # modification positions only in sample 2
both_set = set1 & set2 # modification positions in BOTH samples

# Filter original DataFrames
df1_only = tsv1_df[tsv1_df['Position'].isin(only1_set)]
df2_only = tsv2_df[tsv2_df['Position'].isin(only2_set)]
df_both = tsv1_df[tsv1_df['Position'].isin(both_set)]

df1_only.to_csv(os.path.join('venn_diagram', 'only_'+args.input2+'.tsv'), sep='\t', index=False)
df2_only.to_csv(os.path.join('venn_diagram', 'only_'+args.input4+'.tsv'), sep='\t', index=False)
df_both.to_csv(os.path.join('venn_diagram', 'both_'+args.input2+'n'+arg.input4+'.tsv'), sep='\t', index=False)
```

### Motif analysis and consensus plots
Another method of analysis could be to look at the motifs associated with a given modification. Fortunately, `modkit` includes a program that uses its own machine learning algorithm to tease out different motifs for modifications found in the bedMethyl file. This can easily be done with the `find-motifs` command as follows (note that this is also computationally intensive and therefore should be run on a cluster computer):

```sh
modkit find-motifs -t 12 --in-bedmethyl pileup.bed --ref .fna -o modkit_motifs.tsv
```

This will export a TSV file containing a list of frequently identified motifs for a specific modification. It is worth noting that there are some interesting results where it identifies multiple "variations" of the same motif, however the trailing nucleotide sequences are different and are therefore identified as different. The base notation used by `modkit` also corresponds to the IUPAC nomenclature for bases.

Specifically when investigating methylation modifications, its important to note that these methylation sites typically come in pairs such that there is a methylation site on the forward strand as well as a methylation site on the reverse strand for a given motif. This can be referred to as fully-methylation and hemi-methylation:
* fully-methylation: when there is a methylation site on both the forward and reverse strand for a given motif.
* hemi-methylation: when there is a methylation site on only the forward or reverse strand for a given motif.
The python script below takes the position of each modification position within a modification table TSV file (filtered from the original bedMethyl file) and, using a provided FASTA file, will pull the nucleotide at the given modification position as well as the bases +/- 12 positions away from the modification. This is iterated for each position present in the TSV, save those sequences as a string, and then uses the provided motif to filter  the methylation sites associated to that motif and then identify which of those sites are fully- or hemi-methylated. In order to run this script, the `biopython` python package is required. This can downloaded by the following command:

```sh
pip install biopython
```

The python script for determining the fully- and hemi-methylated sites for a given motif is as follows:

`fully-hemi-methylation.py`
```py
import os
import pandas as pd
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# %% --- Argeparse preamble --- %%

parser = argparse.ArgumentParser(description="Identifies which of the methylation sites are fully-methylated (methylated on both strands) and " \
"hemi-methylated (methylated on one strand) given an annotation table TSV file from the bedMethyl analysis pipeline.")
parser.add_argument("input1", help="Provide the path to the annotation table TSV file, e.g., annotated_table_6mA.tsv")
parser.add_argument("input2", help="Provide the path to the FASTA file (FA or FNA file okay)")
parser.add_argument("input3", help="Provide the motif in a regex format, e.g., 'GATC' or 'CGA[ATCG]{7}TC[CT]'. Note: the motifs MUST be in single quotations!")
parser.add_argument("input4", help="Provide the modification position offset in the motif (in a 0-based indexing system; should be the value given in the modkit motif TSV file)")
args = parser.parse_args()

# Import the TSV input as a DataFrame and only keep the 'Position' and 'Strand' column
tsv_filename = args.input1
consensus_df = pd.read_csv(tsv_filename, sep="\t", usecols=['Position', 'Strand', 'Locus_tag', 'Gene', 'Product'])

# Import the FASTA file path
fna_filename =  args.input2

# Provide motif in regex format
motif = args.input3

# Modification position offset (0-based index)
a_offset = int(args.input4)

# %% --- Preliminary functions to interpret the regex motif --- %%%
# Calculate the length of the motif in regex format
def regex_motif_length(motif):
    length = 0
    i = 0

    while i < len(motif):
        if motif[i] == "[":
            j = motif.index("]", i)
            token_len = 1
            i = j + 1

        elif motif[i] in "ACGT":
            token_len = 1
            i += 1

        else:
            i += 1
            continue

        # check for quantifier
        if i < len(motif) and motif[i] == "{":
            j = motif.index("}", i)
            q = motif[i+1:j]
            n = int(q.split(",")[0])  # {n} or {n,m} → n
            token_len *= n
            i = j + 1

        length += token_len

    return length

# Call the regex_motif_length function to get the length of the motif
motif_len = regex_motif_length(motif)

# Compute reverse complement of a regex motif with [ ] degenerate classes
def reverse_complement_regex(motif):
    comp = str.maketrans("ACGT", "TGCA")
    tokens = []
    i = 0

    while i < len(motif):
        if motif[i] == "[":
            j = motif.index("]", i)
            chars = motif[i+1:j]
            rc_chars = "".join(sorted(c.translate(comp) for c in chars))
            tokens.append(f"[{rc_chars}]")
            i = j + 1

        elif motif[i] == "{":
            j = motif.index("}", i)
            tokens[-1] += motif[i:j+1]   # attach quantifier
            i = j + 1

        else:
            tokens.append(motif[i].translate(comp))
            i += 1

    return "".join(tokens[::-1])


# Call the reverse_complement_regex function to get the reverse compliment strand motif sequence
motif_rc = reverse_complement_regex(motif)

# %% --- Extract the sequence around the each of the modification positions --- %%

# No. bases before and after modification position
window_size = 12 # no. bases before and after modification position

# Load the genome
record = SeqIO.read(fna_filename, "fasta")

# For loop cycling through all the rows in the data frame
for row_index in range(0, consensus_df.shape[0]):

    # Extract the position for that row
    position = consensus_df.loc[row_index, 'Position']

    # Locate the start and end positions for the given window
    start = max(0, position - window_size)
    end = position + window_size - 1

    # If the strand is '-'...
    if consensus_df.loc[row_index, 'Strand'] == '-':
        window_seq = str(record.seq[start:end].reverse_complement())
    else:
        window_seq = str(record.seq[start:end])

    # Insert sequence in the 'Sequence' column of the DataFrame
    consensus_df.loc[row_index, 'Sequence'] = window_seq

# Filter out all the Data Frame such that all of the positions have the same sequence string length
consensus_df = consensus_df[(consensus_df["Sequence"].apply(len) == ((window_size * 2) - 1) )].reset_index(drop=True)

# %% --- Now that we have the sequences for each position, we want to isolate what motifs are associated with each position --- %%

# Calculate the center of the sequence
center = window_size - 1

# Calculate the start and end of the motif using the a_offset
motif_start = center - a_offset
motif_end = motif_start + motif_len

# Filter the consensus_df such the only the positions of the given motif (and reverse complement) are kept
motif_df = consensus_df[
   consensus_df["Sequence"].str.slice(motif_start, motif_end).str.match(motif) |
   consensus_df["Sequence"].str.slice(motif_start, motif_end).str.match(motif_rc)
].reset_index(drop=True)

# Export the Data Frame containing all of the methylation sites associated to the provided motif (if needed)
#motif_df.to_csv('motif_methylation_sites.tsv', sep='\t', index=False)

# Determine and print the number of positions that have the given motif (and reverse complement)
motif_count = motif_df[motif_df['Sequence'].str.slice(motif_start, motif_end).str.match(motif)].shape[0]
motif_rev_count = motif_df[motif_df['Sequence'].str.slice(motif_start, motif_end).str.match(motif_rc)].shape[0]
print('No. of methylation sites: '+ str(consensus_df.shape[0]))
print('No. of sites with the motif ' + motif + ': '+ str(motif_count))
print('No. of sites with the reverse complement (if asymmetrical) motif ' + motif_rc + ': '+ str(motif_rev_count))

# %% --- Now that we have the total number of modification sites associated to a specific motif (and its reverse complement) we want to find the total number of fully- and hemi- methylation sites --- %%

# Counter
i = 0

# Difference in between sister pairs modifications
pos_diff = (motif_len - (a_offset + 1)) - a_offset 

# Empty list for rows to delete
partner_sites = []

# While loop using a counter to index the rows until it reaches the length of the data frame (add the -1 because by proxy we will be looking at the last position)
while i <= motif_df.shape[0]-2:
    # Extract the position for the counter:
    pos1 = motif_df.loc[i, 'Position']
    # Extract the position for the counter + 1 to compare to:
    pos2 = motif_df.loc[i+1, 'Position']
    # Calculate the spacial difference between the two methylation site positions:
    diff = pos2 - pos1

    # If the difference between the positions are 1 or 3 (i.e., GATC or CTSSAG motif) we will record the two positions to remove from the data frame
    if diff == pos_diff:
        # Append the rows that meet this condition (i.e., are "partners") to the list of rows to delete
        partner_sites.append(i)
        partner_sites.append(i+1)

        # Since we have already established that the next position in the data frame is a partner, we can add +2 the counter
        i += 2
    else: 
        i += 1

# Filter the motif data frame with only the fully methylated positions
fully = motif_df.iloc[partner_sites]
print('No. of fully-methylated sites: ' +  str(fully.shape[0]))

# Export the Data Frame containing all of the fully-methylated sites associated to the provided motif (if needed)
#fully.to_csv('fully_methylated_sites.tsv', sep='\t', index=False)
# %% --- Now what if we want to look at the modification sites that are hemi-methylated, i.e., don't have a pair --- %%
hemi = motif_df.drop(partner_sites)
hemi_count = hemi[hemi['Sequence'].str.slice(motif_start, motif_end).str.match(motif)].shape[0]
hemi_rev_count = hemi[hemi['Sequence'].str.slice(motif_start, motif_end).str.match(motif_rc)].shape[0]
print('No. of forward hemi-methylated sites: ' + str(hemi_count))
print('No. of reverse hemi-methylated sites (if asymmetrical): ' + str(hemi_rev_count))

# Export the Data Frame containing all of the hemi-methylated sites associated to the provided motif (if needed)
#hemi.to_csv('hemi_methylated_sites.tsv', sep='\t', index=False)
```
This code will also export the Data Frames as TSV for the methylation sites of the given motif as well as the Data Frames for the fully- and hemi-methylated sites. They are made as comments (denoted by the #) so to have these files exported, simply remove the '#' from the line of code.


While this gives a quantitative way of looking at the motifs present in a strain, a more qualitative way of looking at this would be to generate a consensus logo showing the relative abundance of nucleotide bases present around the modification site. The python script below is an in-house method of generating a consensus logo by using the same methods as `fully-hemi-methylated.py` by taking the sequences surrounding the modification sites. This is iterated for each position present in the TSV, save those sequences as a string, and then uses `logomaker` to generate a consensus plot. In order to run this script, the `logomaker`, as well as the `biopython` and `matplotlib` python packages are required; these can simply be downloaded by the following command:

```py
pip install biopython
pip install logomaker
pip install matplotlib
```

The python script to generate a consensus plot is as follows:

`consensus_logo.py`
```py
import os
import pandas as pd
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logomaker
import matplotlib.pyplot as plt

# %% --- Argeparse preamble --- %%

parser = argparse.ArgumentParser(description="Identifies which of the methylation sites are fully-methylated (methylated on both strands) and " \
"hemi-methylated (methylated on one strand) given an annotation table TSV file from the bedMethyl analysis pipeline.")
parser.add_argument("input1", help="Provide the path to the annotation table TSV file, e.g., annotated_table_6mA.tsv")
parser.add_argument("input2", help="Provide the path to the FASTA file (FA or FNA file okay)")
parser.add_argument("input3", help="Provide the motif in a regex format, e.g., GATC or CGA[ATCG]{7}TC[CT]")
parser.add_argument("input4", help="Provide the modification position offset in the motif (in a 0-based indexing system; should be the value given in the modkit motif TSV file)")
parser.add_argument("input5", help="Provide the strain name, e.g., UA159")
args = parser.parse_args()

# Import the TSV input as a DataFrame and only keep the 'Position' and 'Strand' column
tsv_filename = args.input1
consensus_df = pd.read_csv(tsv_filename, sep="\t", usecols=['Position', 'Strand'])

# Import the FASTA file path
fna_filename =  args.input2

# Provide motif in regex format
motif = args.input3

# Modification position offset (0-based index)
a_offset = int(args.input4)

# Strain name
strain = args.input5

# %% --- Preliminary functions to interpret the regex motif --- %%%
# Calculate the length of the motif in regex format
def regex_motif_length(motif):
    length = 0
    i = 0

    while i < len(motif):
        if motif[i] == "[":
            j = motif.index("]", i)
            token_len = 1
            i = j + 1

        elif motif[i] in "ACGT":
            token_len = 1
            i += 1

        else:
            i += 1
            continue

        # check for quantifier
        if i < len(motif) and motif[i] == "{":
            j = motif.index("}", i)
            q = motif[i+1:j]
            n = int(q.split(",")[0])  # {n} or {n,m} → n
            token_len *= n
            i = j + 1

        length += token_len

    return length

# Call the regex_motif_length function to get the length of the motif
motif_len = regex_motif_length(motif)

# Compute reverse complement of a regex motif with [ ] degenerate classes
def reverse_complement_regex(motif):
    comp = str.maketrans("ACGT", "TGCA")
    tokens = []
    i = 0

    while i < len(motif):
        if motif[i] == "[":
            j = motif.index("]", i)
            chars = motif[i+1:j]
            rc_chars = "".join(sorted(c.translate(comp) for c in chars))
            tokens.append(f"[{rc_chars}]")
            i = j + 1

        elif motif[i] == "{":
            j = motif.index("}", i)
            tokens[-1] += motif[i:j+1]   # attach quantifier
            i = j + 1

        else:
            tokens.append(motif[i].translate(comp))
            i += 1

    return "".join(tokens[::-1])


# Call the reverse_complement_regex function to get the reverse compliment strand motif sequence
motif_rc = reverse_complement_regex(motif)

# %% --- Extract the sequence around the each of the modification positions --- %%

# No. bases before and after modification position
window_size = 12 # no. bases before and after modification position

# Load the genome
record = SeqIO.read(fna_filename, "fasta")

# For loop cycling through all the rows in the data frame
for row_index in range(0, consensus_df.shape[0]):

    # Extract the position for that row
    position = consensus_df.loc[row_index, 'Position']

    # Locate the start and end positions for the given window
    start = max(0, position - window_size)
    end = position + window_size - 1

    # If the strand is '-'...
    if consensus_df.loc[row_index, 'Strand'] == '-':
        window_seq = str(record.seq[start:end].reverse_complement())
    else:
        window_seq = str(record.seq[start:end])

    # Insert sequence in the 'Sequence' column of the DataFrame
    consensus_df.loc[row_index, 'Sequence'] = window_seq

# Filter out all the Data Frame such that all of the positions have the same sequence string length
consensus_df = consensus_df[(consensus_df["Sequence"].apply(len) == ((window_size * 2) - 1) )].reset_index(drop=True)

# %% --- Now that we have the sequences for each position, we want to isolate what motifs are associated with each position --- %%

# Calculate the center of the sequence
center = window_size - 1

# Calculate the start and end of the motif using the a_offset
motif_start = center - a_offset
motif_end = motif_start + motif_len

# Filter the consensus_df such the only the positions of the given motif (and reverse complement) are kept
motif_df = consensus_df[
   consensus_df["Sequence"].str.slice(motif_start, motif_end).str.match(motif) |
   consensus_df["Sequence"].str.slice(motif_start, motif_end).str.match(motif_rc)
].reset_index(drop=True)

# %% --- Create the consensus logo of every methylation site --- %%
sequences = consensus_df['Sequence'].values

# Create a count matrix
counts_df = logomaker.alignment_to_matrix(sequences)

# Plot title
title = "Consensus Logo: " + strain

# Plot logo
plt.figure(figsize=(10, 4))
logo = logomaker.Logo(counts_df)
plt.title(title)
plt.savefig("consensus_logo_" + strain + ".pdf")

# %% --- Create the consensus logo of the motif methylation site --- %%
sequences = motif_df['Sequence'].values

# Create a count matrix
counts_df = logomaker.alignment_to_matrix(sequences)

# Plot title
title = "Motif Consensus Logo: " + strain

# Plot logo
plt.figure(figsize=(10, 4))
logo = logomaker.Logo(counts_df)
plt.title(title)
plt.savefig("motif_consensus_logo_" + strain + ".pdf")
```

### Circos plot generation
Circos is a handy program that can be helpful for generating circular density plots of modification sites. However, it is fairly unintuitive and can be difficult to troubleshoot. Therefore, this section is provided as helpful documentation on how to generate a Circos plot.

1. Download circos to a local computer
    * Go to the [circos website](#https://circos.ca/software/download/) to download the most recent version of circos. For example, at the time of writing this document, the most recent version is `circos-0.69-10.tgz`.
    * Once downloaded, run the following line of code to install it, creating a directory `circos-0.69-10/`:
    ```sh
    tar xzf circos-0.69-10.tgz
    ```
    * Install the Perl dependencies for the appropriate OS:
        - macOS:
        ```sh
        brew install perl cpanminus gd
        ```
    * There may be some Perl modules missing; if needed these can be installed in the following format:
    ```sh
    cpanminus --local-lib=~/perl5 Config::General Font::TTF GD GD::Polyline Clone
    ```
    * Add all the perl dependencies to `~/.zshrc`:
    ```sh
    eval "$(perl -I ~/perl5/lib/perl5 -Mlocal::lib)"
    ```
    * Finally, reload the terminal:
    ```sh
    source ~/.zshrc
    ```
    * Test `circos` to make sure everything runs properly:
    ```sh
    circos-0.69-10/bin/circos -module
    ```
**Note: Everytime you want to run `circos`, you MUST load the shell startup configuration, i.e., `source ~/.zshrc`.**

2. Before running Circos, a few input files are needed to produce the graph:
    * `karyotype.txt`
        - This is a TXT file that is used as the framework for the plot which contains the chromosome name/number, species, and length of the genome. For example, the text below is what the file would look like for *S. mutans*
        
        `chr - smutans smutans 0 2032925 black`

        **Note that this example is of a bacterial genome, i.e., only one chromosome is present and hence the 'chr' name for the chromosome. For a genome containing more than one chromosome (e.g. eukaryotes), numbering the chromosomes is required for the karyotype file to work.**

    * `density.txt`
        - This is a TXT file contains the number of modification sites found within bins of n number of basepairs along the genome. This TXT file is formatted such that it contains a list of `species bin-start bin-end no.-modification-sites`, for example:
        ```
        smutans 0       1000    10
        smutans 1000    2000    9
        smutans 2000    3000    5
        smutans 3000    4000    11
        ...
        ```
        This file can be generated using the python code below that takes a modification table TSV as an input file and automatically generate a `density.txt` file:

`bin_modification.py`
```py
import sys
import csv
from collections import defaultdict

# === USER INPUT ===
input_file = "../../SmuUA159/modifications_tables/filtered_modkit_table_4mC.tsv"   # Your input TSV file
output_file = "./data/4mC_methylation_density.txt"  # Output file for Circos
genome_length = 2032925                     # Replace with your genome length
bin_size = 1000                            # Size of each bin in bp
chromosome_id = "smutans"                  # Must match your karyotype.txt
# ===================

# Prepare bins
num_bins = (genome_length + bin_size - 1) // bin_size
bins = defaultdict(int)

with open(input_file, newline='') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        try:
            pos = int(row[1])  # 2nd column (0-based index 1)
            bin_index = pos // bin_size
            bins[bin_index] += 1
        except (IndexError, ValueError):
            continue  # Skip malformed lines

# Write Circos histogram format
with open(output_file, 'w') as out:
    for i in range(num_bins):
        start = i * bin_size
        end = min(start + bin_size, genome_length)
        count = bins[i]
        out.write(f"{chromosome_id}\t{start}\t{end}\t{count}\n")

print(f"Done! Binned data saved to {output_file}")
```

2. Cont.
    * `circos.conf`
        - Finally, Circos runs and generates a plot via interpreting a configuration file (CONF), which in essensce is a file with text formatted in its own language to generate highly customizable plots. The example below is used to create a density plot of the methylation sites found in *S. mutans*, however this can be used as a template to modify and change parameters according to the specific needs:

```
karyotype = data/karyotype.txt
<<include etc/housekeeping.conf>>
<<include etc/colors_fonts_patterns.conf>>

<image>

image_size = 1600p

</image>

<ideogram>
	# Ideogram position, fill and outline
	radius = 0.75r
	thickness = 40p
	fill = yes
	fill_color = black
	stroke_thickness = 1p
	stroke_color = black
		<spacing>
			default = 0.005r
		</spacing>
		
</ideogram>

show_label = yes
label_font = bold
label_radius   = 0.2r
label_center = yes
label_size = 48

show_ticks = yes
show_tick_labels = yes

<ticks>
	radius = 1.01r
	color = black

		<tick>
			# Major tick marks
			spacing = 200000u
			size = 45p
			thickness = 5p
			show_label = yes
			label_font = bold
			label_size = 40p
			label_offset = 10p
			multiplier = 0.001
			format = %.0f kb
		</tick>

		<tick>
			# Minor tick marks
			spacing = 100000u
			size = 35
			thickness = 5p
		</tick>

</ticks>

<image>
<<include etc/image.conf>>
</image>

<plots>
	<plot>
		type = histogram
		file = data/6mA_methylation_density.txt
		r1 = 0.98r
		r0 = 0.75r
		thickness = 5p
		color = green
		min = 0
		max = 30
			<axes>
				<axis>
					position = 0r
					color = black
					thickness = 5p
				</axis>
				<axis>
					position = 0.25r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.5r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.75r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 1r
					color = black
					thickness = 5p
				</axis>
			</axes>

	</plot>
	<plot>
		type = histogram
		file = data/4mC_methylation_density.txt
		r1 = 0.7r
		r0 = 0.47r
		thickness = 5p
		color = blue
		min = 0
		max = 4
			<axes>
				<axis>
					position = 0r
					color = black
					thickness = 5p
				</axis>
				<axis>
					position = 0.25r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.5r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.75r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 1r
					color = black
					thickness = 5p
				</axis>
			</axes>
	</plot>
	<plot>
		type = histogram
		file = data/5mC_methylation_density.txt
		r1 = 0.42r
		r0 = 0.19r
		thickness = 5p
		color = orange
		min = 0
		max = 4
			<axes>
				<axis>
					position = 0r
					color = black
					thickness = 5p
				</axis>
				<axis>
					position = 0.25r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.5r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 0.75r
					color = lgrey
					thickness = 3p
				</axis>
				<axis>
					position = 1r
					color = black
					thickness = 5p
				</axis>
			</axes>
	</plot>
</plots>
```

3. With the CONF file made, circos can easily be run via the following command:
```sh
circos-0.69-10/bin/circos -conf circos.conf
```
This will export a PNG and a SVG file to the current directory.

## ONT Methylation Pipeline
It is worth noting that the rki-mf1/ONT-methylation pipeline exports a few more files than just the filtered modification tables; it also provides additional files that may be helpful for interpreting the bedMethyl output. For future reference, however, this pipline can be executed by inputting the unmapped BAM file from `dorado` into the following command (assuming all the proper program and packages are downloaded).

### Installing `nextflow`
`nextflow` needs to be downloaded. This can be done in two ways:
 1. Use `curl`: 
 ```sh 
 curl -s https://get.nextflow.io | bash
 ```
 2. Use create a conda environment and use conda to install it:
 ```sh
 mamba create -n ont-methylation_env -c bioconda nextflow
```

With nextflow, the valegale/ONT_methylation analysis pipeline can be downloaded simply by the following:

```sh
nextflow pull valegale/ONT_methylation
# check available release versions and branches:
nextflow info valegale/ONT_methylation
# show the help message for a certain pipeline release version:
nextflow run valegale/ONT_methylation -r 0.0.1 --help
# update the pipeline simply via pulling the code again:
nextflow pull valegale/ONT_methylation
```

### Supplemental downloads for `valegale/ONT_methylation`
Finally, `valegale/ONT_methylation` requires a few additional packages that don't come pre-installed. `valegale/ONT_methylation` requires both `samtools` and `minimap2`. Its recommended to create a new conda environment for these packages (if not done already when downloading `nextflow`).

```sh
conda create -n ont-methylation_env -c bioconda samtools minimap2
```

### Running ont-methylation
**Note that the agrument for the `--fasta` flag must have the same file name as the input `--bam` argument**:

```sh
conda activate ont-methylation_env
./nextflow run valegale/ONT_methylation -r 0.0.1 --fasta barcode00.fna --bam barcode00.bam -profile slurm 
```

The modification tables produced by valgale/ONT-methylation that the [Generating a bedMethyl file using ONT Modkit](#Generating-a-bedMethyl-file-using-ONT-Modkit) section mimics are located in `results/BAM-file-name/modification_tables`.