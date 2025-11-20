#!/usr/bin/env python3
import pandas as pd

#from here: https://github.com/slowkow/allelefrequencies/tree/main

AFND_FILE = r"C:\Users\Kacper\Desktop\immuno\data\afnd.tsv"          # path to the file you just downloaded
OUT_FILE  = r"C:\Users\Kacper\Desktop\immuno\outputs/hla_freq.tsv"      # this is the TSV that select_hla_alleles.py will use

# 1. Load AFND dump
print(f"Loading {AFND_FILE} ...")
df = pd.read_csv(AFND_FILE, sep="\t", dtype=str)

# 2. Keep only HLA genes (group == "hla")
df = df[df["group"] == "hla"].copy()

# 3. Keep only classical loci we care about
#    Class I:  A, B, C
#    Class II: DRB1, DQB1, DPB1
target_genes = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]
df = df[df["gene"].isin(target_genes)].copy()

# 4. Parse allele frequency (alleles_over_2n) as float
df["frequency"] = pd.to_numeric(df["alleles_over_2n"], errors="coerce")
df = df.dropna(subset=["frequency"])

# 5. Build 'locus' column (same as gene)
df["locus"] = df["gene"]

# 6. Ensure allele is in HLA style, e.g. HLA-A*02:01, HLA-DRB1*04:01
def format_allele(row):
    gene = row["gene"]
    allele = row["allele"]
    # Examples in afnd.tsv look like "A*01:01", "DRB1*04:01"
    # We just add the "HLA-" prefix.
    if allele.startswith(gene + "*"):
        return f"HLA-{allele}"
    else:
        # Fallback, just in case the format is slightly different
        return f"HLA-{gene}*{allele}"

df["allele_fmt"] = df.apply(format_allele, axis=1)

# 7. Keep only the columns we need and rename
out = df[["locus", "allele_fmt", "frequency", "population"]].rename(
    columns={"allele_fmt": "allele"}
).copy()

print(df["population"].unique())

# 8. Save to TSV
out.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Saved {len(out)} rows to {OUT_FILE}")
print("Columns:", list(out.columns))
