import pandas as pd
import argparse
import os
from itertools import combinations

def filter_input_data(df):
    """Apply optional filtering to retain significantly changed genes."""
    return df[
        ((df["Log2 fold change"] >= 1) | (df["Log2 fold change"] <= -1)) &
        (df["BY.p.value"] <= 0.05)
    ]

def analyze_gene_presence(input_file, output_dir, apply_filter=False):
    df = pd.read_csv(input_file, sep=",")

    required_cols = ["Gene", "Day"]
    if apply_filter:
        required_cols += ["Log2 fold change", "BY.p.value"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    if apply_filter:
        print("Applying filter: |Log2 fold change| ≥ 1 AND BY.p.value ≤ 0.05")
        df = filter_input_data(df)

    os.makedirs(output_dir, exist_ok=True)

    day_groups = df.groupby("Day")
    day_gene_sets = {}

    for day, group_df in day_groups:
        day_gene_sets[day] = set(group_df["Gene"])

    # Build gene presence matrix
    all_genes = sorted(set.union(*day_gene_sets.values()))
    matrix = pd.DataFrame(index=all_genes)

    for day, genes in day_gene_sets.items():
        matrix[f"Day_{day}"] = [gene in genes for gene in matrix.index]

    matrix = matrix.reset_index().rename(columns={"index": "Gene"})
    presence_file = os.path.join(output_dir, "gene_day_presence_matrix.txt")
    matrix.to_csv(presence_file, sep="\t", index=False)
    print(f"Gene presence matrix saved to: {presence_file}")

    # Build Venn data
    venn_data = []

    # Unique genes per Day
    for day in day_gene_sets:
        other_days = [d for d in day_gene_sets if d != day]
        unique_genes = day_gene_sets[day] - set.union(*(day_gene_sets[d] for d in other_days))
        venn_data.append({
            "Day": str(day),
            "Unique_Genes": ",".join(sorted(unique_genes)),
            "Unique_Count": len(unique_genes)
        })

    # Shared between all Days
    shared_all = set.intersection(*day_gene_sets.values())
    venn_data.append({
        "Day": "Shared_All",
        "Unique_Genes": ",".join(sorted(shared_all)),
        "Unique_Count": len(shared_all)
    })

    # Shared pairwise intersections
    for (day1, day2) in combinations(day_gene_sets.keys(), 2):
        intersection = day_gene_sets[day1] & day_gene_sets[day2]
        if intersection:  # Only report non-empty intersections
            venn_data.append({
                "Day": f"{day1}_{day2}",
                "Unique_Genes": ",".join(sorted(intersection)),
                "Unique_Count": len(intersection)
            })

    venn_file = os.path.join(output_dir, "gene_venn_data.txt")
    pd.DataFrame(venn_data).to_csv(venn_file, sep="\t", index=False)
    print(f"Venn diagram gene data saved to: {venn_file}")

def main():
    parser = argparse.ArgumentParser(description="Analyze gene presence by Day for Venn diagram generation.")
    parser.add_argument("input_file", help="Path to input tab-delimited file with 'Gene' and 'Day' columns.")
    parser.add_argument("--output_dir", "-o", default="output_day_venn", help="Directory to store output files.")
    parser.add_argument("--filter", action="store_true", help="Apply filtering: |Log2 fold change| >= 1 and BY.p.value <= 0.05")

    args = parser.parse_args()

    try:
        analyze_gene_presence(args.input_file, args.output_dir, args.filter)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
