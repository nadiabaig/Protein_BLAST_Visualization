#!/usr/bin/env python3
"""
Protein BLAST (-outfmt 6) → plots + tables (robust, chunked Excel).

Example:
  python3 prot_blast.py --in results.tsv --outdir out --topn 5 --pdf --xlsx

Recommended BLAST outfmt (for richer plots):
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids sscinames salltitles"
"""

import argparse
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

STANDARD_12 = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]
EXCEL_MAX_ROWS = 1_048_576
CHUNK_ROWS = 1_000_000         
MAX_SHEETNAME = 31


def ensure_dir(p): os.makedirs(p, exist_ok=True)

def sanitize_filename(s: str) -> str:
    bad = '<>:"/\\|?* \t\n\r'
    out = "".join(ch if ch not in bad else "_" for ch in str(s))
    return out[:120]

def save_plot(name: str, outdir: str, pdf: bool):
    f_png = os.path.join(outdir, f"{name}.png")
    f_pdf = os.path.join(outdir, f"{name}.pdf")
    plt.tight_layout()
    plt.savefig(f_png, dpi=300, bbox_inches="tight")
    if pdf:
        plt.savefig(f_pdf, bbox_inches="tight")
    plt.close()
    return [f_png] + ([f_pdf] if pdf else [])


def read_blast_outfmt6(path: str) -> pd.DataFrame:
    # Try headered
    try:
        df = pd.read_csv(path, sep="\t", engine="python")
        if {"qseqid","sseqid"}.issubset(df.columns):
            return df
    except Exception:
        pass
    # Fallback: headerless
    df = pd.read_csv(path, sep="\t", header=None, engine="python")
    n = df.shape[1]
    if n < 12:
        raise ValueError(f"{path}: expected ≥12 columns for -outfmt 6, got {n}.")
    if n == 12:
        df.columns = STANDARD_12
    else:
        df.columns = STANDARD_12 + [f"extra{i}" for i in range(1, n-12+1)]
    return df

def coerce_numeric(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    coerce_numeric(df, ["pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send",
                        "evalue","bitscore","qlen","slen"])
    # coverage
    if "qlen" in df.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["qcov"] = (df["length"] / df["qlen"]) * 100
    else:
        df["qcov"] = np.nan
    # log10 evalue
    if "evalue" in df.columns:
        df["log10_evalue"] = np.log10(df["evalue"].replace(0, np.finfo(float).tiny))
    else:
        df["log10_evalue"] = np.nan
    return df

# ---------- tables ----------
def top_hits_per_query(df: pd.DataFrame, topn: int) -> pd.DataFrame:
    keys, asc = [], []
    if "evalue" in df.columns: keys.append("evalue"); asc.append(True)
    if "bitscore" in df.columns: keys.append("bitscore"); asc.append(False)
    if keys:
        s = df.sort_values(by=["qseqid"]+keys, ascending=[True]+asc)
    else:
        s = df.sort_values(by=["qseqid"])
    return s.groupby("qseqid", as_index=False).head(topn)

def per_query_summary(df: pd.DataFrame) -> pd.DataFrame:
    grp = df.groupby("qseqid")
    def best(sub):
        s = sub.sort_values(by=["evalue","bitscore"], ascending=[True, False]).iloc[0]
        return pd.Series({
            "best_sseqid": s.get("sseqid", np.nan),
            "best_sscinames": s.get("sscinames", np.nan),
            "best_title": s.get("salltitles", np.nan),
            "best_pident": s.get("pident", np.nan),
            "best_length": s.get("length", np.nan),
            "best_evalue": s.get("evalue", np.nan),
            "best_bitscore": s.get("bitscore", np.nan),
            "best_qcov": s.get("qcov", np.nan)
        })
    summ = grp.apply(best).reset_index()
    extra = grp.agg(n_hits=("sseqid","count"),
                    median_pident=("pident","median"),
                    median_qcov=("qcov","median")).reset_index()
    return summ.merge(extra, on="qseqid", how="left")

def per_taxon_summary(df: pd.DataFrame) -> pd.DataFrame:
    if "sscinames" not in df.columns: return pd.DataFrame()
    g = df.groupby("sscinames")
    return g.agg(n_hits=("sseqid","count"),
                 best_identity=("pident","max"),
                 best_coverage=("qcov","max"),
                 best_bitscore=("bitscore","max")).reset_index().sort_values("n_hits", ascending=False)

# ---------- plots ----------
def plot_hist_pident(df, outdir, pdf):
    if "pident" not in df.columns: return []
    plt.figure()
    plt.hist(df["pident"].dropna(), bins=30)
    plt.xlabel("% identity"); plt.ylabel("Count"); plt.title("Distribution of % identity")
    return save_plot("hist_pident", outdir, pdf)

def plot_hist_qcov(df, outdir, pdf):
    if "qcov" not in df.columns or df["qcov"].dropna().empty: return []
    plt.figure()
    plt.hist(df["qcov"].dropna(), bins=30)
    plt.xlabel("Query coverage (%)"); plt.ylabel("Count"); plt.title("Distribution of query coverage")
    return save_plot("hist_qcov", outdir, pdf)

def plot_scatter_pident_len(df, outdir, pdf):
    if "pident" not in df.columns or "length" not in df.columns: return []
    x = df["length"].astype(float); y = df["pident"].astype(float)
    m = (~x.isna()) & (~y.isna())
    if not m.any(): return []
    plt.figure()
    plt.scatter(x[m], y[m], s=8, alpha=0.6)
    plt.xlabel("Alignment length (aa)"); plt.ylabel("% identity"); plt.title("% identity vs alignment length")
    return save_plot("scatter_pident_length", outdir, pdf)

def plot_bar_taxa(df, outdir, pdf, top=20):
    if "sscinames" not in df.columns: return []
    counts = df["sscinames"].value_counts().head(top)
    if counts.empty: return []
    plt.figure()
    counts.plot(kind="bar")
    plt.ylabel("Count of hits"); plt.title(f"Taxonomic distribution (top {top})")
    return save_plot("bar_taxa", outdir, pdf)

def plot_coverage_map(df, outdir, pdf, query=None):
    need = {"qlen","qstart","qend","qseqid"}
    if not need.issubset(df.columns): return []
    if query is None:
        query = df["qseqid"].value_counts().index[0]
    sub = df[df["qseqid"] == query].copy()
    if sub.empty: return []
    qlen = pd.to_numeric(sub["qlen"], errors="coerce").dropna()
    if qlen.empty: return []
    qlen = int(qlen.iloc[0])

    plt.figure(figsize=(8, 2.2))
    plt.hlines(y=0.5, xmin=0, xmax=qlen, linewidth=3)
    for _, r in sub.iterrows():
        try:
            qs = float(r["qstart"]); qe = float(r["qend"])
            x0, x1 = sorted([qs, qe])
            if not (math.isnan(x0) or math.isnan(x1)):
                plt.hlines(y=0.5, xmin=x0, xmax=x1, linewidth=1)
        except Exception:
            continue
    plt.xlabel(f"Query position (aa; 0..{qlen})")
    plt.yticks([]); plt.title(f"Coverage map: {query}")
    return save_plot(f"coverage_{sanitize_filename(query)}", outdir, pdf)


def _clean_for_excel(df: pd.DataFrame) -> pd.DataFrame:
    def _fix(x):
        if isinstance(x, str):
            x = "".join(ch if (ord(ch) >= 32 or ch in ("\n","\t")) else " " for ch in x)
            if len(x) > 32000:
                x = x[:32000]
        return x
    return df.applymap(_fix)

def _safe_sheet(base: str, idx=None):
    name = base if idx is None else f"{base}_{idx}"
    return name[:MAX_SHEETNAME]

def _write_df_chunked(writer, df: pd.DataFrame, base_sheet: str):
    df = _clean_for_excel(df)
    n = len(df)
    if n <= CHUNK_ROWS:
        df.to_excel(writer, index=False, sheet_name=_safe_sheet(base_sheet))
        return
    start, part = 0, 1
    while start < n:
        stop = min(start + CHUNK_ROWS, n)
        df.iloc[start:stop].to_excel(writer, index=False,
                                     sheet_name=_safe_sheet(base_sheet, part))
        start = stop
        part += 1

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Protein BLAST visualization & tables (-outfmt 6).")
    ap.add_argument("--in", dest="infile", required=True, help="BLAST tabular file (-outfmt 6).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--topn", type=int, default=5, help="Top N hits per query.")
    ap.add_argument("--query", type=str, default=None, help="Query ID for coverage map (optional).")
    ap.add_argument("--pdf", action="store_true", help="Also write vector PDF plots.")
    ap.add_argument("--xlsx", action="store_true", help="Also write Excel workbook (with chunking).")
    args = ap.parse_args()

    ensure_dir(args.outdir)

    df = read_blast_outfmt6(args.infile)
    df = add_derived(df)

    # --- tables (CSV) ---
    all_hits_csv = os.path.join(args.outdir, "all_hits_clean.csv")
    df.to_csv(all_hits_csv, index=False)

    top_hits = top_hits_per_query(df, args.topn)
    top_csv = os.path.join(args.outdir, "top_hits_per_query.csv")
    top_hits.to_csv(top_csv, index=False)

    q_summary = per_query_summary(df)
    qs_csv = os.path.join(args.outdir, "query_summary.csv")
    q_summary.to_csv(qs_csv, index=False)

    tax_summary = per_taxon_summary(df)
    if not tax_summary.empty:
        ts_csv = os.path.join(args.outdir, "taxon_summary.csv")
        tax_summary.to_csv(ts_csv, index=False)

    # --- Excel (optional, chunked) ---
    if args.xlsx:
        try:
            engine = "xlsxwriter"
            try:
                import xlsxwriter  # noqa
            except Exception:
                engine = "openpyxl"
            xlsx_path = os.path.join(args.outdir, "results.xlsx")
            with pd.ExcelWriter(xlsx_path, engine=engine) as xl:
                _write_df_chunked(xl, df, "all_hits")
                _write_df_chunked(xl, top_hits, "top_hits")
                _write_df_chunked(xl, q_summary, "query_summary")
                if not tax_summary.empty:
                    _write_df_chunked(xl, tax_summary, "taxon_summary")
            print(f"Excel workbook written: {xlsx_path}")
        except Exception as e:
            print(f"[warn] Excel write failed: {e}\n       (CSV files were written successfully.)")

    # --- plots ---
    written = []
    written += plot_hist_pident(df, args.outdir, args.pdf)
    written += plot_hist_qcov(df, args.outdir, args.pdf)
    written += plot_scatter_pident_len(df, args.outdir, args.pdf)
    written += plot_bar_taxa(df, args.outdir, args.pdf)
    written += plot_coverage_map(df, args.outdir, args.pdf, query=args.query)

    # --- report ---
    print("CSV tables:")
    print(" -", all_hits_csv)
    print(" -", top_csv)
    print(" -", qs_csv)
    if not tax_summary.empty:
        print(" -", ts_csv)
    print("\nPlots:")
    if written:
        for p in written: print(" -", p)
    else:
        print(" - (no plots generated; check columns)")

if __name__ == "__main__":
    main()
