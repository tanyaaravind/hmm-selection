"""
Data preparation for 1000 Genomes ABO/Rh loci.

This script focuses on Wellington's responsibilities:
- Read VCF (or pre-extracted allele frequency table).
- Compute per-population allele frequencies and DeltaAF = |AF_popA - AF_popB|.
- Save tidy tables for downstream HMM/benchmarking.

The code uses cyvcf2 if available; otherwise it can operate on a small
CSV allele-frequency export. Keep inputs small when testing locally.
"""

import os
import sys
import json
import pandas as pd
from typing import Dict, List, Tuple

try:
    from cyvcf2 import VCF  # Optional, not in base requirements
except ImportError:
    VCF = None


def load_sample_map(sample_map_path: str) -> Dict[str, str]:
    """
    Load mapping of sample ID -> population code (e.g., YRI, CEU, CHB).
    Expected CSV with columns: sample_id,pop.
    """
    df = pd.read_csv(sample_map_path)
    if not {"sample_id", "pop"}.issubset(df.columns):
        raise ValueError("sample_map must have columns sample_id,pop")
    return dict(zip(df["sample_id"], df["pop"]))


def allele_frequencies_from_vcf(vcf_path: str, sample_map: Dict[str, str]) -> pd.DataFrame:
    """
    Compute allele frequencies per population from a VCF file.

    Returns a tidy DataFrame with columns:
    chrom, pos, ref, alt, pop, alt_count, total_alleles, af
    """
    if VCF is None:
        raise ImportError("cyvcf2 is required for VCF parsing. Install with `pip install cyvcf2`.")

    vcf = VCF(vcf_path)
    sample_to_pop = sample_map

    # Build population -> indices lookup
    pop_to_indices: Dict[str, List[int]] = {}
    for idx, sample in enumerate(vcf.samples):
        pop = sample_to_pop.get(sample)
        if pop is None:
            continue
        pop_to_indices.setdefault(pop, []).append(idx)

    records = []
    for variant in vcf:
        for pop, idxs in pop_to_indices.items():
            # GTs are tuples like (0,1)
            alt_count = 0
            total = 0
            for sample_idx in idxs:
                gt = variant.genotypes[sample_idx]
                # genotype is (allele1, allele2, phased flag)
                if len(gt) < 2:
                    continue
                a1, a2 = gt[0], gt[1]
                if a1 >= 0:
                    total += 1
                    alt_count += a1
                if a2 >= 0:
                    total += 1
                    alt_count += a2
            if total == 0:
                continue
            af = alt_count / total
            records.append(
                {
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": variant.ALT[0],
                    "pop": pop,
                    "alt_count": alt_count,
                    "total_alleles": total,
                    "af": af,
                }
            )

    return pd.DataFrame.from_records(records)


def delta_af_table(af_df: pd.DataFrame, pop_a: str, pop_b: str) -> pd.DataFrame:
    """
    Compute DeltaAF = |AF_pop_a - AF_pop_b| for matched variants.
    """
    a = af_df[af_df["pop"] == pop_a][["chrom", "pos", "af"]].rename(columns={"af": f"af_{pop_a}"})
    b = af_df[af_df["pop"] == pop_b][["chrom", "pos", "af"]].rename(columns={"af": f"af_{pop_b}"})
    merged = a.merge(b, on=["chrom", "pos"], how="inner")
    merged["delta_af"] = (merged[f"af_{pop_a}"] - merged[f"af_{pop_b}"]).abs()
    return merged


def save_tables(af_df: pd.DataFrame, delta_df: pd.DataFrame, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    af_path = os.path.join(out_dir, "allele_frequencies.csv")
    delta_path = os.path.join(out_dir, "delta_af.csv")
    af_df.to_csv(af_path, index=False)
    delta_df.to_csv(delta_path, index=False)
    return af_path, delta_path


def main():
    """
    Minimal CLI for quick tests.

    Examples:
    python src/data_prep.py --vcf data/abo_slice.vcf.gz --sample-map data/samples.csv --pop-a YRI --pop-b CEU --out results/abo_freqs

    For environments without VCF parsing, you can provide a precomputed
    allele frequency CSV via --af-table.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Compute allele frequencies and DeltaAF for 1000G loci.")
    parser.add_argument("--vcf", type=str, help="Path to VCF (bgzipped) for ABO/Rh slice.")
    parser.add_argument("--sample-map", type=str, help="CSV with sample_id,pop mapping.")
    parser.add_argument("--af-table", type=str, help="Optional: existing allele frequency CSV (chrom,pos,pop,af).")
    parser.add_argument("--pop-a", type=str, default="YRI")
    parser.add_argument("--pop-b", type=str, default="CEU")
    parser.add_argument("--out", type=str, default="results/abo_freqs")
    args = parser.parse_args()

    if args.af_table:
        af_df = pd.read_csv(args.af_table)
    else:
        if not args.vcf or not args.sample_map:
            parser.error("Provide either --af-table or both --vcf and --sample-map")
        sample_map = load_sample_map(args.sample_map)
        af_df = allele_frequencies_from_vcf(args.vcf, sample_map)

    delta_df = delta_af_table(af_df, args.pop_a, args.pop_b)
    af_path, delta_path = save_tables(af_df, delta_df, args.out)
    print(json.dumps({"allele_frequencies": af_path, "delta_af": delta_path}, indent=2))


if __name__ == "__main__":
    main()
