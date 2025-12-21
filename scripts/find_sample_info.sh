#!/bin/bash

echo "Searching for sample metadata files on 1000 Genomes FTP..."
echo ""

PATHS=(
    "/vol1/ftp/technical/working/20130606_sample_info"
    "/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    "/vol1/ftp/phase3/data"
    "/vol1/ftp/technical/working"
)

for path in "${PATHS[@]}"; do
    echo "Trying: ftp://ftp.1000genomes.ebi.ac.uk$path"
    curl -s "ftp://ftp.1000genomes.ebi.ac.uk$path/" 2>/dev/null | grep -i "ped\|panel\|sample" | head -3
    echo ""
done

echo "Alternative: Try listing the technical directory:"
curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/" 2>/dev/null | head -10

