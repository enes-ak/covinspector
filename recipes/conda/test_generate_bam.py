"""Generate a minimal synthetic BAM for conda build tests."""
import os
import pysam

header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr17", "LN": 83257441}],
    "RG": [{"ID": "test", "SM": "TEST_SAMPLE", "PL": "ILLUMINA"}],
}

bam_path = "test_synthetic.bam"
with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
    for i in range(200):
        a = pysam.AlignedSegment()
        a.query_name = f"read_{i}"
        a.query_sequence = "A" * 150
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 43044294 + (i * 400)
        a.mapping_quality = 60
        a.cigar = [(0, 150)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 150)
        a.set_tag("RG", "test")
        outf.write(a)

pysam.sort("-o", "test_sorted.bam", bam_path)
pysam.index("test_sorted.bam")
os.remove(bam_path)
print("Synthetic BAM created.")
