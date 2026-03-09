import pandas as pd
from collections import defaultdict

# List of input files to process
input_files = {
    "Bp": "genes_unique_Bp.txt",  
    "Yp": "genes_unique_Yp.txt",  
    "Shared": "genes_shared_between_bugs.txt"
}

# Dictionary to hold counts: {gene_set: {group: count}}
gene_set_counts = defaultdict(lambda: defaultdict(int))

for group, filename in input_files.items():
    try:
        df = pd.read_csv(filename, sep="\t")
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        continue

    # Drop rows with missing Gene.sets values
    df = df.dropna(subset=["Gene.sets"])

    for sets in df["Gene.sets"]:
        # Split by comma and strip whitespace
        for gene_set in map(str.strip, sets.split(",")):
            if gene_set:  # skip empty
                gene_set_counts[gene_set][group] += 1

# Convert to DataFrame
summary_df = pd.DataFrame.from_dict(gene_set_counts, orient="index")
summary_df.index.name = "Gene.set"
summary_df = summary_df.fillna(0).astype(float).reset_index()

# Divide the 'Shared' column by (number of input files - 1)
num_inputs_minus_one = len(input_files) - 1
if "Shared" in summary_df.columns and num_inputs_minus_one > 0:
    summary_df["Shared"] = summary_df["Shared"] / num_inputs_minus_one

# Convert all numeric columns back to int (rounding first)
numeric_cols = summary_df.columns.drop("Gene.set")
summary_df[numeric_cols] = summary_df[numeric_cols].round().astype(int)

# Sort by Gene.set alphabetically
summary_df = summary_df.sort_values(by="Gene.set")

# Output to file
summary_df.to_csv("gene_set_summary_counts.txt", sep="\t", index=False)

print("Summary file written to 'gene_set_summary_counts.txt'")
