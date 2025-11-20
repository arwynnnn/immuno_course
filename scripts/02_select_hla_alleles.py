#!/usr/bin/env python
"""
Select representative HLA alleles for NetMHCpan (Class I) and NetMHCIIpan (Class II)
based on allele frequency tables (e.g. derived from AFND), restricted to alleles
actually supported by the NetMHC tools.

Input
-----
A TSV/CSV file with at least the following columns:
  - locus     : e.g. A, B, C, DRB1, DQB1, DPB1
  - allele    : full allele name, e.g. HLA-A*02:01, HLA-DRB1*04:01
  - frequency : allele frequency (0–1) in the population of interest

Optional columns:
  - population : if present, you can filter by population name via --population.

Additionally, provide the NetMHC allele list files:
  - NetMHCpan-4.1:   MHC_allele_names.txt      (one allele per line)
  - NetMHCIIpan-4.3: alleles_name.txt         (tabular, DR/DQ/DP columns)

Usage
-----
python select_hla_alleles.py \
    --freq-table hla_freq.tsv \
    --population "Northern Europe" \
    --max-classI 20 \
    --max-classII 15 \
    --coverage 0.9 \
    --netmhcpan-alleles MHC_allele_names.txt \
    --netmhciipan-alleles alleles_name.txt

The script prints two tables to stdout and also writes:
  - selected_classI_alleles.txt  (NetMHCpan format, e.g. HLA-A02:01)
  - selected_classII_alleles.txt (NetMHCIIpan format, e.g. DRB1*04:01, DPB1*04:01)
"""

import argparse
import sys
import pandas as pd


CLASS_I_LOCI = ["A", "B", "C"]
CLASS_II_LOCI = ["DRB1", "DQB1", "DPB1"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Select representative HLA alleles for NetMHCpan/NetMHCIIpan"
    )
    p.add_argument(
        "--freq-table",
        required=True,
        help="TSV/CSV file with columns: locus, allele, frequency [, population]",
    )
    p.add_argument(
        "--sep", default="\t", help="Field separator for freq table (default: TAB)"
    )
    p.add_argument(
        "--population",
        nargs="+",
        default=None,
        help="Population label(s) to filter on; use multiple values to combine",
    )
    p.add_argument(
        "--coverage",
        type=float,
        default=0.9,
        help="Target cumulative frequency per locus (default: 0.9)",
    )
    p.add_argument(
        "--max-classI",
        type=int,
        default=20,
        help="Maximum number of class I alleles to return (default: 20)",
    )
    p.add_argument(
        "--max-classII",
        type=int,
        default=15,
        help="Maximum number of class II alleles to return (default: 15)",
    )
    p.add_argument(
        "--netmhcpan-alleles",
        required=True,
        help="NetMHCpan allele list file (MHC_allele_names.txt)",
    )
    p.add_argument(
        "--netmhciipan-alleles",
        required=True,
        help="NetMHCIIpan allele list file (alleles_name.txt)",
    )
    return p.parse_args()


def select_for_class(df, loci, max_alleles, coverage):
    """
    Select alleles for a given class (I or II).

    Strategy:
      - For each locus in `loci`:
        - Collapse duplicates per allele (keep highest frequency)
        - Sort alleles by descending frequency
        - Add alleles until per-locus cumulative frequency >= `coverage`
      - If total selected > max_alleles, keep the top `max_alleles` by frequency.

    Returns a DataFrame with columns at least: locus, allele, frequency
    """
    selected_rows = []

    for locus in loci:
        sub = df[df["locus"] == locus].copy()
        if sub.empty:
            continue

        # make sure frequency is numeric
        sub["frequency"] = sub["frequency"].astype(float)

        # ✅ collapse multiple rows of the same allele (keep highest freq)
        sub = sub.sort_values("frequency", ascending=False)
        sub = sub.drop_duplicates(subset=["allele"], keep="first")

        cum = 0.0
        for _, row in sub.iterrows():
            selected_rows.append(row)
            cum += float(row["frequency"])
            if cum >= coverage:
                break

    if not selected_rows:
        return pd.DataFrame(columns=df.columns)

    sel_df = pd.DataFrame(selected_rows)

    # If we exceeded max_alleles, keep the top by frequency
    sel_df["frequency"] = sel_df["frequency"].astype(float)
    if sel_df.shape[0] > max_alleles:
        sel_df = sel_df.sort_values("frequency", ascending=False).head(max_alleles)

    return sel_df.reset_index(drop=True)


# ---------- Formatting helpers ----------


def format_netmhcpan(allele):
    """
    Convert AFND-style class I allele (e.g. HLA-A*02:01 or A*02:01)
    to NetMHCpan style used in MHC_allele_names.txt, e.g. HLA-A02:01.
    """
    a = str(allele).strip()
    if not a:
        return a

    # Ensure "HLA-" prefix is present
    if not a.startswith("HLA-"):
        # e.g. "A*02:01" -> "HLA-A*02:01"
        if a[0] in {"A", "B", "C"} and (a[1] == "*" or a[1] == "0"):
            a = "HLA-" + a
        else:
            # Unexpected format; return unchanged
            return a

    # Remove the '*' between locus and numbers
    a = a.replace("*", "")
    return a


def format_netmhciipan(allele):
    """
    Convert AFND-style class II allele (e.g. HLA-DRB1*04:01)
    into NetMHCIIpan allele name as used in alleles_name.txt, e.g. DRB1*04:01.
    (Simply drop the 'HLA-' prefix if present; keep ':' and '*')
    """
    a = str(allele).strip()
    if a.startswith("HLA-"):
        a = a[4:]
    return a


# ---------- Load NetMHC allele lists ----------


def load_netmhcpan_supported(path):
    """
    Load NetMHCpan allele names from MHC_allele_names.txt and
    return a set of supported human class I alleles (HLA-A/B/C...).
    """
    supported = set()
    with open(path, "r") as f:
        for line in f:
            name = line.strip()
            if not name:
                continue
            if name.startswith("HLA-A") or name.startswith("HLA-B") or name.startswith(
                "HLA-C"
            ):
                supported.add(name)
    return supported


def load_netmhciipan_supported(path):
    """
    Load NetMHCIIpan allele names from alleles_name.txt and return
    three sets: supported DRB1, DQB1, and DPB1 alleles
    (naming as in the file, e.g. DRB1*04:01, DQB1*03:01, DPB1*04:01).
    """
    drb1 = set()
    dqb1 = set()
    dpb1 = set()

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("DR "):
                continue
            parts = line.split()
            for tok in parts:
                if tok.startswith("DRB1*"):
                    drb1.add(tok)
                elif tok.startswith("DQB1*"):
                    dqb1.add(tok)
                elif tok.startswith("DPB1*"):
                    dpb1.add(tok)

    return drb1, dqb1, dpb1

def format_netmhciipan_output(allele):
    """
    Convert e.g. HLA-DRB1*07:01 or DRB1*07:01 to DRB1_0701
    for NetMHCIIpan text input.
    """
    a = str(allele).strip()
    if not a:
        return a

    # drop HLA- prefix if present
    if a.startswith("HLA-"):
        a = a[4:]

    if "*" in a:
        locus, rest = a.split("*", 1)
    else:
        # fall back: no *, just remove colons and replace remaining * by _
        return a.replace(":", "").replace("*", "_")

    # remove colon between fields, e.g. 07:01 -> 0701
    rest = rest.replace(":", "")
    return f"{locus}_{rest}"


def main():
    args = parse_args()

    # Load NetMHC-supported alleles
    supported_I = load_netmhcpan_supported(args.netmhcpan_alleles)
    supported_DRB1, supported_DQB1, supported_DPB1 = load_netmhciipan_supported(
        args.netmhciipan_alleles
    )

    df = pd.read_csv(args.freq_table, sep=args.sep, dtype=str)

    required_cols = {"locus", "allele", "frequency"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(
            f"ERROR: freq table must contain columns: {', '.join(sorted(required_cols))}"
        )

    # cast frequency to float
    df["frequency"] = df["frequency"].astype(float)

    # Filter populations
    if args.population is not None:
        if "population" not in df.columns:
            sys.exit(
                "ERROR: --population was given but 'population' column is missing in freq table"
            )

        df = df[df["population"].isin(args.population)].copy()

        if df.empty:
            sys.exit(f"ERROR: no rows matching populations: {args.population}")

    # ----- Class I -----
    classI_df = df[df["locus"].isin(CLASS_I_LOCI)].copy()
    # Map to NetMHCpan naming and keep only supported alleles
    if not classI_df.empty:
        classI_df["netmhcpan_allele"] = classI_df["allele"].map(format_netmhcpan)
        before = classI_df.shape[0]
        classI_df = classI_df[classI_df["netmhcpan_allele"].isin(supported_I)].copy()
        after = classI_df.shape[0]
        if after == 0:
            print(
                "WARNING: No Class I alleles from freq table are supported by NetMHCpan with given allele list.",
                file=sys.stderr,
            )
    else:
        classI_df["netmhcpan_allele"] = []

    sel_I = select_for_class(classI_df, CLASS_I_LOCI, args.max_classI, args.coverage)
    # Ensure mapping exists on the selected subset
    if not sel_I.empty and "netmhcpan_allele" not in sel_I.columns:
        sel_I["netmhcpan_allele"] = sel_I["allele"].map(format_netmhcpan)

    # ----- Class II -----
    classII_df = df[df["locus"].isin(CLASS_II_LOCI)].copy()

    if not classII_df.empty:
        classII_df["netmhciipan_allele"] = classII_df["allele"].map(format_netmhciipan)

        def is_supported(row):
            allele = row["netmhciipan_allele"]
            locus = row["locus"]
            if locus == "DRB1":
                return allele in supported_DRB1
            elif locus == "DQB1":
                return allele in supported_DQB1
            elif locus == "DPB1":
                return allele in supported_DPB1
            else:
                return False

        before_II = classII_df.shape[0]
        classII_df = classII_df[classII_df.apply(is_supported, axis=1)].copy()
        after_II = classII_df.shape[0]
        if after_II == 0:
            print(
                "WARNING: No Class II alleles from freq table are supported by NetMHCIIpan with given allele list.",
                file=sys.stderr,
            )
    else:
        classII_df["netmhciipan_allele"] = []

    sel_II = select_for_class(classII_df, CLASS_II_LOCI, args.max_classII, args.coverage)
    if not sel_II.empty and "netmhciipan_allele" not in sel_II.columns:
        sel_II["netmhciipan_allele"] = sel_II["allele"].map(format_netmhciipan)

    # ---------- Write outputs (comma-separated, no spaces) ----------
    classI_file = "selected_classI_alleles.txt"
    classII_file = "selected_classII_alleles.txt"

    with open(classI_file, "w") as f:
        if sel_I.empty:
            f.write("")
        else:
            f.write(",".join(sel_I["netmhcpan_allele"].astype(str)))

    with open(classII_file, "w") as f:
        formatted = [format_netmhciipan_output(a) for a in sel_II["netmhciipan_allele"]]
        f.write(",".join(formatted))

    # ---------- Console summary ----------
    print(f"Selected Class I alleles (n={sel_I.shape[0]}):")
    if not sel_I.empty:
        print(sel_I[["locus", "allele", "frequency", "netmhcpan_allele"]])
    else:
        print("(none)")
    print()

    print(f"Selected Class II alleles (n={sel_II.shape[0]}):")
    if not sel_II.empty:
        print(sel_II[["locus", "allele", "frequency", "netmhciipan_allele"]])
    else:
        print("(none)")
    print()

    print(f"Wrote Class I allele list to {classI_file}")
    print(f"Wrote Class II allele list to {classII_file}")


if __name__ == "__main__":
    main()
