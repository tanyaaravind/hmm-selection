
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Iterable, Tuple


def _safe_array(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError("Expected 1D array-like input")
    return arr


def compute_pairwise_fst(
    alt_counts_a: Iterable[float],
    total_alleles_a: Iterable[float],
    alt_counts_b: Iterable[float],
    total_alleles_b: Iterable[float],
) -> np.ndarray:
    
    alt_a = _safe_array(alt_counts_a)
    alt_b = _safe_array(alt_counts_b)
    tot_a = _safe_array(total_alleles_a)
    tot_b = _safe_array(total_alleles_b)

    if not (len(alt_a) == len(alt_b) == len(tot_a) == len(tot_b)):
        raise ValueError("Input arrays must be the same length.")

    p1 = np.divide(alt_a, tot_a, out=np.zeros_like(alt_a, dtype=float), where=tot_a > 0)
    p2 = np.divide(alt_b, tot_b, out=np.zeros_like(alt_b, dtype=float), where=tot_b > 0)
    p_bar = (p1 * tot_a + p2 * tot_b) / np.maximum(tot_a + tot_b, 1e-9)

    num = (p1 - p2) ** 2 - (p1 * (1 - p1) / np.maximum(tot_a - 1, 1e-9)) - (
        p2 * (1 - p2) / np.maximum(tot_b - 1, 1e-9)
    )
    den = p_bar * (1 - p_bar)

    fst = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    fst[fst < 0] = 0.0 
    return fst


def fst_from_af_table(af_df: pd.DataFrame, pop_a: str, pop_b: str) -> pd.DataFrame:

    required_cols = {"chrom", "pos", "pop", "af"}
    if not required_cols.issubset(af_df.columns):
        missing = ", ".join(sorted(required_cols - set(af_df.columns)))
        raise ValueError(f"af_df is missing columns: {missing}")

    subset = af_df[af_df["pop"].isin([pop_a, pop_b])]
    pivot = subset.pivot_table(index=["chrom", "pos"], columns="pop", values="af")
    pivot = pivot.dropna(subset=[pop_a, pop_b])
    pivot = pivot.reset_index()

    has_counts = {"alt_count", "total_alleles"}.issubset(af_df.columns)
    if has_counts:
        count_subset = af_df[af_df["pop"].isin([pop_a, pop_b])]
        count_pivot_alt = count_subset.pivot_table(index=["chrom", "pos"], columns="pop", values="alt_count")
        count_pivot_total = count_subset.pivot_table(index=["chrom", "pos"], columns="pop", values="total_alleles")
        count_pivot_alt = count_pivot_alt.reindex(pivot.set_index(["chrom", "pos"]).index)
        count_pivot_total = count_pivot_total.reindex(pivot.set_index(["chrom", "pos"]).index)

        fst_vals = compute_pairwise_fst(
            count_pivot_alt[pop_a].values,
            count_pivot_total[pop_a].values,
            count_pivot_alt[pop_b].values,
            count_pivot_total[pop_b].values,
        )
    else:
        p1 = pivot[pop_a].values
        p2 = pivot[pop_b].values
        p_bar = 0.5 * (p1 + p2)
        denom = 2 * p_bar * (1 - p_bar)
        fst_vals = np.divide((p1 - p2) ** 2, denom, out=np.zeros_like(p1), where=denom > 0)

    pivot["fst"] = fst_vals
    return pivot[["chrom", "pos", "fst"]]


def tajimas_d_from_counts(alt_counts: Iterable[int], total_alleles: Iterable[int]) -> float:

    alt = _safe_array(alt_counts)
    total = _safe_array(total_alleles)

    if len(alt) != len(total):
        raise ValueError("alt_counts and total_alleles must be the same length.")

    n = np.round(np.mean(total)).astype(int)
    if n < 2:
        raise ValueError("Need at least two chromosomes to compute Tajima's D.")

    freqs = np.divide(alt, total, out=np.zeros_like(alt, dtype=float), where=total > 0)
    segregating = (freqs > 0) & (freqs < 1)
    s = segregating.sum()
    if s == 0:
        return 0.0

    pi = np.sum(2 * freqs[segregating] * (1 - freqs[segregating]))

    a1 = np.sum(1 / np.arange(1, n))
    a2 = np.sum(1 / (np.arange(1, n) ** 2))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    denom = np.sqrt(e1 * s + e2 * s * (s - 1))
    if denom == 0:
        return 0.0
    return (pi - s / a1) / denom


def ihs_from_ehh(
    positions: Iterable[float],
    ehh_ref: Iterable[float],
    ehh_alt: Iterable[float],
) -> Tuple[float, float, float]:

    pos = _safe_array(positions)
    ref = _safe_array(ehh_ref)
    alt = _safe_array(ehh_alt)

    if not (len(pos) == len(ref) == len(alt)):
        raise ValueError("positions, ehh_ref, ehh_alt must be the same length.")
    if not np.all(np.diff(pos) >= 0):
        raise ValueError("positions must be sorted.")

    i_ref = np.trapz(ref, pos)
    i_alt = np.trapz(alt, pos)
    if i_ref <= 0 or i_alt <= 0:
        return i_ref, i_alt, 0.0

    ihs = np.log(i_ref / i_alt)
    return i_ref, i_alt, ihs


__all__ = [
    "compute_pairwise_fst",
    "fst_from_af_table",
    "tajimas_d_from_counts",
    "ihs_from_ehh",
]
