import pandas as pd
from collections import defaultdict
import argparse
import os

def filter_input_data(df):
    """Apply optional filtering to retain significantly changed genes."""
    return df[
        ((df["Log2 fold change"] >= 1) | (df["Log2 fold change"] <= -1)) &
        (df["BY.p.value"] <= 0.05)
    ]

def split_by_bug_and_gene(input_file, output_dir, apply_filter):
    df = pd.read_csv(input_file, sep="\t")

    required_cols = ["Gene", "Bug", "Log2 fold change", "BY.p.value"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    if apply_filter:
        print("🔍 Applying filtering: |Log2 fold change| >= 1 AND BY.p.value <= 0.05")
        df = filter_input_data(df)

    bug_groups = df.groupby("Bug")
    bug_gene_sets = {}

    for bug, group_df in bug_groups:
        bug_gene_sets[bug] = set(group_df["Gene"])
        group_df.to_csv(os.path.join(output_dir, f"genes_unique_{bug}.txt"), sep="\t", index=False)

    if len(bug_gene_sets) == 2:
        bug1, bug2 = list(bug_gene_sets.keys())
        shared_genes = bug_gene_sets[bug1] & bug_gene_sets[bug2]
        unique_genes_bug1 = bug_gene_sets[bug1] - shared_genes
        unique_genes_bug2 = bug_gene_sets[bug2] - shared_genes

        df[(df["Bug"] == bug1) & (df["Gene"].isin(unique_genes_bug1))].to_csv(
            os.path.join(output_dir, f"genes_unique_{bug1}.txt"), sep="\t", index=False)
        df[(df["Bug"] == bug2) & (df["Gene"].isin(unique_genes_bug2))].to_csv(
            os.path.join(output_dir, f"genes_unique_{bug2}.txt"), sep="\t", index=False)
        df[df["Gene"].isin(shared_genes)].to_csv(
            os.path.join(output_dir, "genes_shared_between_bugs.txt"), sep="\t", index=False)

        return bug1, bug2
    else:
        raise ValueError("Exactly two Bug groups are required.")

def count_gene_sets(bug1, bug2, output_dir):
    input_files = {
        bug1: os.path.join(output_dir, f"genes_unique_{bug1}.txt"),
        bug2: os.path.join(output_dir, f"genes_unique_{bug2}.txt"),
        "Shared": os.path.join(output_dir, "genes_shared_between_bugs.txt")
    }

    gene_set_counts = defaultdict(lambda: defaultdict(int))

    for group, filename in input_files.items():
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found. Skipping.")
            continue

        df = pd.read_csv(filename, sep="\t")
        if "Gene.sets" not in df.columns:
            print(f"Warning: File '{filename}' missing 'Gene.sets' column. Skipping.")
            continue

        df = df.dropna(subset=["Gene.sets"])

        for sets in df["Gene.sets"]:
            for gene_set in map(str.strip, sets.split(",")):
                if gene_set:
                    gene_set_counts[gene_set][group] += 1

    summary_df = pd.DataFrame.from_dict(gene_set_counts, orient="index")
    summary_df.index.name = "Gene.set"
    summary_df = summary_df.fillna(0).astype(int).reset_index()
    summary_df = summary_df.sort_values(by="Gene.set")

    summary_file = os.path.join(output_dir, "gene_set_summary_counts.txt")
    summary_df.to_csv(summary_file, sep="\t", index=False)
    print(f"✅ Output written to {summary_file}")

def main():
    parser = argparse.ArgumentParser(description="Split gene data by bug and summarize gene sets.")
    parser.add_argument("input_file", help="Path to the tab-delimited input file.")
    parser.add_argument("--output_dir", "-o", default="output", help="Directory to save output files (default: ./output)")
    parser.add_argument("--filter", action="store_true", help="Apply filtering: |Log2 fold change| >= 1 and BY.p.value <= 0.05")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        bug1, bug2 = split_by_bug_and_gene(args.input_file, args.output_dir, args.filter)
        count_gene_sets(bug1, bug2, args.output_dir)
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
