#!/bin/bash

cd "$(dirname "$0")/../data/raw_vcf" || exit 1

echo "Attempting to download chromosome 9 VCF from 1000 Genomes..."
echo ""

echo "Trying v5b version..."
if wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz; then
    echo "✅ Successfully downloaded v5b VCF file!"
    VERSION="v5b"
else
    echo "v5b failed, trying v5a..."
    if wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; then
        echo "✅ Successfully downloaded v5a VCF file!"
        VERSION="v5a"
    else
        echo "❌ Both versions failed. Listing available files..."
        curl -s ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ | grep "chr9"
        exit 1
    fi
fi

echo ""
echo "Downloading index file..."
if wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_${VERSION}.20130502.genotypes.vcf.gz.tbi; then
    echo "✅ Successfully downloaded index file!"
else
    echo "⚠️  Index file download failed, but VCF file is available"
fi

echo ""
echo "Verifying download..."
if [ -f "ALL.chr9.phase3_shapeit2_mvncall_integrated_${VERSION}.20130502.genotypes.vcf.gz" ]; then
    FILE_SIZE=$(ls -lh "ALL.chr9.phase3_shapeit2_mvncall_integrated_${VERSION}.20130502.genotypes.vcf.gz" | awk '{print $5}')
    echo "✅ File exists: ALL.chr9.phase3_shapeit2_mvncall_integrated_${VERSION}.20130502.genotypes.vcf.gz"
    echo "   Size: $FILE_SIZE"
    echo ""
    echo "Next steps:"
    echo "  1. Extract ABO region: bcftools view -r 9:133240000-133290000 ALL.chr9...vcf.gz -Oz -o ../../abo_region.vcf.gz"
    echo "  2. Process with: python src/data_prep.py --vcf data/abo_region.vcf.gz --sample-map data/samples.csv --pop-a YRI --pop-b CEU --out results/abo_expanded"
else
    echo "❌ File not found after download attempt"
    exit 1
fi

