import sys
import csv
from collections import defaultdict

# === USER INPUT ===
input_file = "../CGAN7TCY_methylation_sites.tsv"   # Your input TSV file
output_file = "./data/UA159_CGAN7TCY_methylation_density.txt"  # Output file for Circos
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
            pos = int(row[0])  # 2nd column (0-based index 1)
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

print(f"✅ Done! Binned data saved to {output_file}")
