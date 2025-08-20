import pandas as pd
from collections import defaultdict

# List of input files to process
input_files = {
    "Bug1": "genes_unique_Bp.txt",  
    "Bug2": "genes_unique_Yp.txt",  
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
summary_df = summary_df.fillna(0).astype(int).reset_index()

# Sort by Gene.set alphabetically
summary_df = summary_df.sort_values(by="Gene.set")

# Output to file
summary_df.to_csv("gene_set_summary_counts.txt", sep="\t", index=False)

print("✅ Summary file written to 'gene_set_summary_counts.txt'")
