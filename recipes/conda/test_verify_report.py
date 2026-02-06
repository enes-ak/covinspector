"""Verify report contains expected sections."""
with open("test_report.md") as f:
    content = f.read()
    assert "# covsnap Coverage Report" in content
    assert "BRCA1" in content
    assert "GENCODE v44" in content or "gencode_v44" in content
    assert "hg38" in content
    print("All report checks passed.")
