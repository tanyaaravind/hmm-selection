#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd
import requests
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def annotate_positions(chrom: str, positions: pd.Series, pause=0.2):
    rows = []
    for pos in positions:
        region = f"{chrom}:{int(pos)}-{int(pos)}"
        url = f"{ENSEMBL_REST}/overlap/region/human/{region}?feature=gene"
        try:
            r = requests.get(url, headers=HEADERS, timeout=10)
            if r.status_code == 200:
                data = r.json()
                if len(data) == 0:
                    rows.append({"pos": int(pos), "annotated": False, "gene_id": None, "gene_name": None, "biotype": None, "start": None, "end": None, "strand": None, "annotation_status": "no-gene"})
                else:
                    g = data[0]
                    rows.append({
                        "pos": int(pos),
                        "annotated": True,
                        "gene_id": g.get("id"),
                        "gene_name": g.get("external_name"),
                        "biotype": g.get("biotype"),
                        "start": g.get("start"),
                        "end": g.get("end"),
                        "strand": g.get("strand"),
                        "annotation_status": "ok",
                    })
            else:
                rows.append({"pos": int(pos), "annotated": False, "annotation_status": f"http-{r.status_code}"})
        except Exception as e:
            rows.append({"pos": int(pos), "annotated": False, "annotation_status": "api-failed"})
        time.sleep(pause)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fst", required=True, help="Path to per-site FST CSV (columns: chrom,pos,fst)")
    parser.add_argument("--top", type=int, default=20, help="Top N sites to annotate by FST")
    parser.add_argument("--chrom", default="9", help="Chromosome name (no chr prefix needed)")
    parser.add_argument("--out", default=None, help="Output CSV path (optional)")
    args = parser.parse_args()

    fst = pd.read_csv(args.fst)
    top = fst.sort_values("fst", ascending=False).head(args.top)
    pos_series = top["pos"].astype(int)

    print(f"Annotating top {len(pos_series)} positions using Ensembl REST... (chrom {args.chrom})")
    ann = annotate_positions(args.chrom, pos_series)

    out_path = Path(args.out) if args.out else OUT / (Path(args.fst).stem + f"_top{args.top}_annotated.csv")
    merged = pd.merge(top.reset_index(drop=True), ann, on="pos", how="left")
    merged.to_csv(out_path, index=False)
    print(f"Wrote annotations to {out_path}")


if __name__ == "__main__":
    main()
