"""Verify raw TSV has required columns."""
with open("test_raw.tsv") as f:
    lines = []
    for l in f:
        if not l.startswith("#covsnap"):
            lines.append(l)

header = lines[0].strip().split("\t")
required = (
    "target_id", "contig", "start", "end", "length_bp",
    "mean_depth", "median_depth", "min_depth", "max_depth",
    "pct_zero", "pct_ge_20", "engine_used", "build",
)
for col in required:
    assert col in header, f"Missing column: {col}"
assert len(lines) >= 2, "Expected at least header + 1 data row"
data = lines[1].strip().split("\t")
assert data[header.index("target_id")] == "BRCA1"
assert data[header.index("build")] == "hg38"
print("All TSV checks passed.")
