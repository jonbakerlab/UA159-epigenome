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
parser.add_argument("input3", help="Provide the motif in a regex format, e.g., 'GATC' or 'CGANNNNNNTCY'. Note: the motifs MUST be in single quotations!")
parser.add_argument("input4", help="Provide the modification position offset in the motif (in a 0-based indexing system; should be the value given in the modkit motif TSV file)")
args = parser.parse_args()

# Import the TSV input as a DataFrame and only keep the 'Position' and 'Strand' column
tsv_filename = args.input1
df = pd.read_csv(tsv_filename, sep="\t", usecols=['Contig', 'Position', 'Strand', 'Total_coverage', 'Percent_modified', 'Locus_tag', 'Gene', 'Product'])

# Import the FASTA file path
fna_filename =  args.input2

# Provide motif in regex format
motif_iupac = args.input3

# Modification position offset (0-based index)
a_offset = int(args.input4)

# %% --- Preliminary functions to interpret the regex motif --- %%%

IUPAC_TO_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"
}

def iupac_to_regex(motif):
    return "".join(IUPAC_TO_REGEX[base] for base in motif.upper())

motif = iupac_to_regex(motif_iupac)

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
records = SeqIO.to_dict(SeqIO.parse(fna_filename, "fasta"))

# For loop cycling through all the rows in the data frame
for row_index in range(0, df.shape[0]):

    # Extract the position for that row
    position = df.loc[row_index, 'Position']
    strand = df.loc[row_index, 'Strand']
    contig = df.loc[row_index, 'Contig']

    record = records[contig]

    # Locate the start and end positions for the given window
    start = max(0, position - window_size)
    end = position + window_size - 1

    # If the strand is '-'...
    if strand == '-':
        window_seq = str(record.seq[start:end].reverse_complement())
    else:
        window_seq = str(record.seq[start:end])

    # Insert sequence in the 'Sequence' column of the DataFrame
    df.loc[row_index, 'Sequence'] = window_seq

# %% --- Now that we have the sequences for each position, we want to isolate what the positions that belong to the motif of interest --- %%
# Create a list of both the forward and reverse (complement) motifs and format them into a regex format via re.compile
motif_patterns = [
    re.compile(motif),
    re.compile(motif_rc)
]

# The methylation site is ALWAYS in the center of the sequence extracted, i.e., position 11
CENTER = 11

# Create a function to identify if the motif (fwd & rc) surrounds the methylated A
def has_centered_motif(seq):
    for pattern in motif_patterns:
        if any(m.start() <= CENTER < m.end()
               for m in pattern.finditer(seq)):
            return True
    return False

# Apply the function to filter out all methylation sites that do not belong to the motif of interest
df_motif = df[df['Sequence'].apply(has_centered_motif)].copy()

# %% -- This section is simply to count how many of the remaing sites associated to the motif of interest belong to the fwd or rc motif -- %%

# Create a function (similar to the one above) that takes a SINGLE motif and determines if it is centered around the A
def has_centered_match(seq, pattern):
    return any(
        m.start() <= CENTER < m.end()
        for m in pattern.finditer(seq)
    )

# Apply the function inputting the fwd motif
forward_mask = df_motif["Sequence"].apply(
    lambda seq: has_centered_match(seq, re.compile(motif))
)

# Apply the finction inputting the rc motif
reverse_mask = df_motif["Sequence"].apply(
    lambda seq: has_centered_match(seq, re.compile(motif_rc))
)

# Determine the count (number of) fwd- and rc-associated sites
forward_count = forward_mask.sum()
reverse_count = (reverse_mask & ~forward_mask).sum()
total_count = forward_count + reverse_count

# Determine and print the number of positions that have the given motif (and reverse complement)
print('No. of methylation sites: '+ str(df.shape[0]))
print(f'No. of sites with the motif: {total_count}')
print(f'Forward motif ({motif}): {forward_count}')
print(f'Reverse complement motif ({motif_rc}): {reverse_count}')

# %% --- Now we want to determine the number of fully- and hemi- methylated sites --- %%
# First split the df_motif into forward and reverse strand data frames
df_fwd = df_motif[df_motif['Strand'] == '+'].copy()
df_rev = df_motif[df_motif['Strand'] == '-'].copy()

# We want to calculate the difference that would be found between two fully-methylated sites
pos_diff = abs( (motif_len - (a_offset + 1)) - a_offset )

# Take each position in df_fwd and df_rev, +/- the theoretical fully-methylated sites position difference, and search the opposite strand positions to see if that position exists.
# If it does, it is fully-methylated and if it doesn't, it is hemi-methylated.

# Looking at the forward strand
rev_positions = set(df_rev["Position"])

def find_partner_fwd(pos):
    if pos + pos_diff in rev_positions:
        return pos + pos_diff
    elif pos - pos_diff in rev_positions:
        return pos - pos_diff
    return None

df_fwd["partner_position"] = df_fwd["Position"].apply(find_partner_fwd)
df_fwd["fully_methylated"] = df_fwd["partner_position"].notna()

# Looking at the reverse strand
fwd_positions = set(df_fwd["Position"])

def find_partner_rev(pos):
    if pos + pos_diff in fwd_positions:
        return pos + pos_diff
    elif pos - pos_diff in fwd_positions:
        return pos - pos_diff
    return None

df_rev["partner_position"] = df_rev["Position"].apply(find_partner_rev)
df_rev["fully_methylated"] = df_rev["partner_position"].notna()

# Tallying up the totals from the results
fwd_fully = (df_fwd["fully_methylated"]).sum()
rev_fully = (df_rev["fully_methylated"]).sum()
fwd_hemi = (~df_fwd["fully_methylated"]).sum()
rev_hemi = (~df_rev["fully_methylated"]).sum()
total_fully = fwd_fully + rev_fully
percent_fully = (total_fully/total_count) * 100
total_hemi = fwd_hemi + rev_hemi
percent_hemi = (total_hemi/total_count) * 100


print(f'No. of fully-methylated sites: {total_fully} ({percent_fully: .1f}%)')
print(f'No. of hemi-methylated sites : {total_hemi} ({percent_hemi: .1f}%)')

# %% --- With those values made, now we just want recombine and reorder the df_fwd and df_rev so that we can separate them now into dfs containing either only fully- or hemi- methylated sites --- %%
df_combined = pd.concat([df_fwd, df_rev]).sort_index()
df_fully = df_combined[df_combined["fully_methylated"]].copy()
df_hemi = df_combined[~df_combined["fully_methylated"]].copy()

# %% Exporting Data Frames:

# Export the Data Frame containing all of the methylation sites and their sequences (if needed)
#df.to_csv('methylation_sites_sequences.tsv', sep='\t', index=False)

# Export the Data Frame containing all of the methylation sites associated  to the motif of interest
#df_motif.to_csv('motif_methylation_sites.tsv', sep='\t', index=False)

# Export the Data Frames containing all the fully- or hemi-methylated sites associated to the motif of interest
#df_fully.to_csv('fully-methylated_sites.tsv', sep='\t', index=False)
#df_hemi.to_csv('hemi-methylated_sites.tsv', sep='\t', index=False)