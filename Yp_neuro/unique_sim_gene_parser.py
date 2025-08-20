import pandas as pd

# Load the tab-delimited file
input_file = "/Users/briansmith/Downloads/test_file.txt"  # Replace with your actual filename
df = pd.read_csv(input_file, sep="\t")

# Ensure "Gene" and "Bug" columns exist
assert "Gene" in df.columns and "Bug" in df.columns, "Required columns missing."

# Split data by Bug
bug_groups = df.groupby("Bug")

# Dictionary to hold gene sets for comparison
bug_gene_sets = {}

# Save separate files for each Bug group
for bug, group_df in bug_groups:
    bug_gene_sets[bug] = set(group_df["Gene"])
    output_file = f"genes_unique_{bug}.txt"
    # We'll temporarily write all genes, will update this below after filtering
    group_df.to_csv(output_file, sep="\t", index=False)

# Find shared and unique genes
if len(bug_gene_sets) == 2:
    bug1, bug2 = bug_gene_sets.keys()
    shared_genes = bug_gene_sets[bug1] & bug_gene_sets[bug2]
    unique_genes_bug1 = bug_gene_sets[bug1] - shared_genes
    unique_genes_bug2 = bug_gene_sets[bug2] - shared_genes

    # Filter and write unique genes for each Bug
    df[(df["Bug"] == bug1) & (df["Gene"].isin(unique_genes_bug1))].to_csv(f"genes_unique_{bug1}.txt", sep="\t", index=False)
    df[(df["Bug"] == bug2) & (df["Gene"].isin(unique_genes_bug2))].to_csv(f"genes_unique_{bug2}.txt", sep="\t", index=False)

    # Save shared genes to a separate file
    df[df["Gene"].isin(shared_genes)].to_csv("genes_shared_between_bugs.txt", sep="\t", index=False)
else:
    print("Warning: More than two Bugs detected. Shared gene comparison is only done between two Bugs.")

print("Done!")