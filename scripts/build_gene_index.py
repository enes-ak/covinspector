#!/usr/bin/env python3
"""Build the hg38 gene index from GENCODE GTF.

This script is run ONCE at package build time (not at runtime).
It produces:
  - data/hg38_genes.tsv.gz      (bgzipped, tabix-indexable)
  - data/hg38_genes.tsv.gz.tbi  (tabix index)
  - data/hg38_gene_aliases.json.gz (case-folded alias map)

Usage:
    python scripts/build_gene_index.py \
        --gtf gencode.v44.annotation.gtf.gz \
        --out-dir data/

Prerequisites:
    - Download gencode.v44.annotation.gtf.gz from GENCODE (one-time).
    - Requires: pysam (for bgzip/tabix), or tabix CLI.

The GTF is NOT shipped with the package. Only the derived index files are.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
from collections import defaultdict


def parse_gtf_genes(gtf_path):
    """Extract gene-level records from GENCODE GTF.

    Yields dicts with: contig, start, end, strand, gene_name, gene_id, gene_type.
    Coordinates are converted to 0-based half-open.
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            contig = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based → 0-based
            end = int(fields[4])         # GTF end is 1-based inclusive → 0-based half-open is same value
            strand = fields[6]

            # Parse attributes
            attrs = {}
            for attr in fields[8].rstrip(";").split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, _, val = attr.partition(" ")
                attrs[key] = val.strip('"')

            gene_name = attrs.get("gene_name", "")
            gene_id = attrs.get("gene_id", "")
            gene_type = attrs.get("gene_type", "")

            if not gene_name:
                continue

            yield {
                "contig": contig,
                "start": start,
                "end": end,
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "gene_type": gene_type,
            }


def parse_gtf_exons(gtf_path):
    """Extract exon-level records for MANE Select / canonical transcripts.

    Yields dicts with: contig, start, end, exon_id, gene_name, exon_number, transcript_id.
    """
    # First pass: identify MANE_Select transcripts (or fall back to longest)
    # For simplicity, extract all exons tagged with tag "MANE_Select"
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "exon":
                continue

            attrs_raw = fields[8]
            # Quick check for MANE_Select tag
            if "MANE_Select" not in attrs_raw:
                continue

            contig = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])

            attrs = {}
            for attr in attrs_raw.rstrip(";").split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, _, val = attr.partition(" ")
                attrs[key] = val.strip('"')

            yield {
                "contig": contig,
                "start": start,
                "end": end,
                "exon_id": attrs.get("exon_id", ""),
                "gene_name": attrs.get("gene_name", ""),
                "exon_number": attrs.get("exon_number", "0"),
                "transcript_id": attrs.get("transcript_id", ""),
            }


def build_alias_map(genes):
    """Build a case-folded alias → canonical gene_name map."""
    alias_map = {}
    for g in genes:
        name = g["gene_name"]
        # Map lowercase to canonical
        alias_map[name.lower()] = name
        # Could add HGNC alias lookups here in future
    return alias_map


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, help="GENCODE GTF (gzipped)")
    parser.add_argument("--out-dir", default="data/", help="Output directory")
    parser.add_argument("--build-exons", action="store_true",
                        help="Also build exon index")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ---- Gene index ----
    print(f"Parsing genes from {args.gtf}...", file=sys.stderr)
    genes = list(parse_gtf_genes(args.gtf))
    print(f"  Found {len(genes)} gene records.", file=sys.stderr)

    # Sort by contig (chr-sorted) and start position
    chrom_order = {f"chr{i}": i for i in range(1, 23)}
    chrom_order["chrX"] = 23
    chrom_order["chrY"] = 24
    chrom_order["chrM"] = 25
    genes.sort(key=lambda g: (chrom_order.get(g["contig"], 99), g["start"]))

    # Write uncompressed TSV, then bgzip + tabix
    tsv_path = os.path.join(args.out_dir, "hg38_genes.tsv")
    gz_path = tsv_path + ".gz"

    with open(tsv_path, "w") as f:
        f.write("#contig\tstart\tend\tgene_name\tgene_id\tstrand\tgene_type\n")
        for g in genes:
            f.write(f"{g['contig']}\t{g['start']}\t{g['end']}\t"
                    f"{g['gene_name']}\t{g['gene_id']}\t{g['strand']}\t"
                    f"{g['gene_type']}\n")

    # bgzip
    subprocess.run(["bgzip", "-f", tsv_path], check=True)
    print(f"  Written: {gz_path}", file=sys.stderr)

    # tabix
    subprocess.run(["tabix", "-s1", "-b2", "-e3", gz_path], check=True)
    print(f"  Indexed: {gz_path}.tbi", file=sys.stderr)

    # ---- Alias map ----
    alias_map = build_alias_map(genes)
    alias_path = os.path.join(args.out_dir, "hg38_gene_aliases.json.gz")
    with gzip.open(alias_path, "wt") as f:
        json.dump(alias_map, f, separators=(",", ":"))
    print(f"  Written: {alias_path} ({len(alias_map)} entries)", file=sys.stderr)

    # ---- Exon index (optional) ----
    if args.build_exons:
        print(f"Parsing exons (MANE_Select) from {args.gtf}...", file=sys.stderr)
        exons = list(parse_gtf_exons(args.gtf))
        print(f"  Found {len(exons)} exon records.", file=sys.stderr)

        exons.sort(key=lambda e: (chrom_order.get(e["contig"], 99), e["start"]))

        exon_tsv = os.path.join(args.out_dir, "hg38_exons.bed")
        exon_gz = exon_tsv + ".gz"

        with open(exon_tsv, "w") as f:
            f.write("#contig\tstart\tend\texon_id\tgene_name\texon_number\ttranscript_id\n")
            for e in exons:
                f.write(f"{e['contig']}\t{e['start']}\t{e['end']}\t"
                        f"{e['exon_id']}\t{e['gene_name']}\t{e['exon_number']}\t"
                        f"{e['transcript_id']}\n")

        subprocess.run(["bgzip", "-f", exon_tsv], check=True)
        subprocess.run(["tabix", "-s1", "-b2", "-e3", exon_gz], check=True)
        print(f"  Written: {exon_gz} + .tbi", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
