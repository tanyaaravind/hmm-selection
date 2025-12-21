#!/bin/bash

echo "Searching for chromosome 9 VCF files on 1000 Genomes FTP..."
echo ""

echo "Available chromosome 9 files:"
curl -s ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ | \
  grep -i "chr9" | \
  awk '{print $NF}' | \
  sort

echo ""
echo "Most likely filename:"
curl -s ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ | \
  grep -i "chr9.*vcf.gz" | \
  head -1 | \
  awk '{print "  " $NF}'

echo ""
echo "To download, use:"
echo "  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/<filename_from_above>"

