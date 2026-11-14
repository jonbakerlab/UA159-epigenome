import os
import pandas as pd
import argparse

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Filter the bedMethyl file exported by modkit and create individual TSV files for each modification")
parser.add_argument("input1", help="Provide the path to the bedMethyl files, e.g. pileup.bed")
args = parser.parse_args()

# %% Filtering thresholds
total_coverage = 10
percent_modified = 0.50

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
# Create a dictionary for the methylation modifications and their assiciated IDs
# Add more modifications if needed
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