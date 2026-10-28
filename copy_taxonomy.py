import argparse

import pandas as pd

parser = argparse.ArgumentParser(
    description="Copy phage taxonomy from database genomes onto the other genomes in the same MOP-UP subgroup"
)
parser.add_argument("input", help="MOP-UP Master.csv file (e.g. nameofrunMaster.csv)")
parser.add_argument("output", help="Path for the taxonomy-filled csv file")
args = parser.parse_args()

df = pd.read_csv(args.input)

taxonomy_cols = [
    "phage_Class",
    "phage_Order",
    "phage_Family",
    "MOP-UP Cluster",
    "phage_Genus"
]

for subgroup, group in df.groupby("Subgroup"):

    for col in taxonomy_cols:

        nonempty = group[col].dropna().astype(str)
        nonempty = nonempty[nonempty.str.strip() != ""]

        unique_vals = nonempty.unique()

        if len(unique_vals) > 1:
            print(
                f"Warning: subgroup {subgroup} has multiple "
                f"{col} values: {list(unique_vals)}"
            )

        if len(unique_vals) > 0:
            fill_value = unique_vals[0]

            mask = (
                (df["Subgroup"] == subgroup)
                & (
                    df[col].isna()
                    | (df[col].astype(str).str.strip() == "")
                )
            )

            df.loc[mask, col] = fill_value

df.to_csv(args.output, index=False)