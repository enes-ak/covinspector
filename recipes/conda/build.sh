#!/bin/bash
set -euo pipefail

# Install the Python package
${PYTHON} -m pip install . --no-deps --no-build-isolation -vvv

# Ensure packaged data files are in the right place
# The gene/exon index files are included in the sdist/wheel via
# package_data in pyproject.toml / setup.cfg.
# Verify they exist after install:
DATA_DIR=$(${PYTHON} -c "import covsnap; import os; print(os.path.join(os.path.dirname(covsnap.__file__), 'data'))")
echo "Data directory: ${DATA_DIR}"
if [ ! -f "${DATA_DIR}/hg38_genes.tsv.gz" ]; then
    echo "WARNING: hg38_genes.tsv.gz not found — full gene index not available (built-in fallback will be used)" >&2
fi

echo "covsnap installed successfully."
