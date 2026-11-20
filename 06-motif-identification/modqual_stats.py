#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():

    parser = argparse.ArgumentParser(
        description="Calculate statistics and plot a histogram of mod_qual values from modkit extract full output."
    )

    parser.add_argument(
        "extract_file",
        help="Output file from 'modkit extract full'"
    )

    args = parser.parse_args()

    # Read file
    df = pd.read_csv(args.extract_file, sep="\t")

    # Keep only 6mA calls
    if "mod_code" in df.columns:
        df = df[df["mod_code"] == "a"].copy()

    if len(df) == 0:
        print("No 6mA calls found.")
        return

    mod_qual = df["mod_qual"].dropna()

    # Statistics
    mean_val = mod_qual.mean()
    median_val = mod_qual.median()
    std_val = mod_qual.std()

    print("\nmod_qual Statistics")
    print("-" * 40)
    print(f"N Calls : {len(mod_qual):,}")
    print(f"Mean    : {mean_val:.2f}")
    print(f"Median  : {median_val:.2f}")
    print(f"Std Dev : {std_val:.2f}")
    print(f"Min     : {mod_qual.min():.2f}")
    print(f"Max     : {mod_qual.max():.2f}")

    # Histogram
    plt.figure(figsize=(8, 6))

    plt.hist(
        mod_qual,
        bins=50,
        edgecolor="black",
        alpha=0.75
    )

    plt.axvline(
        mean_val,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean = {mean_val:.2f}"
    )

    plt.axvline(
        median_val,
        color="blue",
        linestyle="--",
        linewidth=2,
        label=f"Median = {median_val:.2f}"
    )

    plt.xlabel("Modification Quality Score (mod_qual)")
    plt.ylabel("Count")
    plt.title("Distribution of mod_qual Scores")
    plt.legend()

    plt.tight_layout()

    output_file = "mod_qual_histogram.png"
    plt.savefig(output_file, dpi=300)

    print(f"\nHistogram written to: {output_file}")


if __name__ == "__main__":
    main()
