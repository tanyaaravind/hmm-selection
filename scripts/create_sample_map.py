#!/usr/bin/env python3
"""
Convert 1000 Genomes PED file to sample map CSV format for data_prep.py

Usage:
    python scripts/create_sample_map.py
"""

import pandas as pd
import os
import sys

# Get project root
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)

def create_sample_map(ped_path=None, output_path=None):
    """
    Convert 1000 Genomes PED file to sample map CSV.
    
    Parameters:
    -----------
    ped_path : str
        Path to 20130606_g1k.ped file
    output_path : str
        Path to save samples.csv
    """
    if ped_path is None:
        ped_path = os.path.join(project_root, 'data', '20130606_g1k.ped')
    if output_path is None:
        output_path = os.path.join(project_root, 'data', 'samples.csv')
    
    # Try panel file first (more reliable)
    panel_path = os.path.join(project_root, 'data', 'integrated_call_samples_v3.20130502.ALL.panel')
    if os.path.exists(panel_path):
        print(f"Using panel file: {panel_path}")
        panel = pd.read_csv(panel_path, sep='\t')
        sample_map = panel[['sample', 'pop']].copy()
        sample_map.columns = ['sample_id', 'pop']
        sample_map.to_csv(output_path, index=False)
        
        print(f"\n✅ Created sample map: {output_path}")
        print(f"   Total samples: {len(sample_map)}")
        print(f"   Populations: {len(sample_map['pop'].unique())}")
        print(f"\n   Population breakdown:")
        pop_counts = sample_map['pop'].value_counts().sort_index()
        for pop, count in pop_counts.items():
            print(f"     {pop}: {count} samples")
        return output_path
    
    if not os.path.exists(ped_path):
        print(f"Error: Neither panel nor PED file found!")
        print("\nDownload one of these files:")
        print("  Option 1 (recommended):")
        print("    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
        print("  Option 2:")
        print("    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
        sys.exit(1)
    
    print(f"Reading PED file: {ped_path}")
    
    # Read the PED file (tab-separated)
    try:
        ped = pd.read_csv(
            ped_path, 
            sep='\t', 
            header=None,
            names=['FamilyID', 'SampleID', 'PaternalID', 'MaternalID', 
                   'Sex', 'Phenotype', 'Population']
        )
    except Exception as e:
        print(f"Error reading PED file: {e}")
        print("Make sure the file is tab-separated and has the correct format")
        sys.exit(1)
    
    # Extract just SampleID and Population
    sample_map = ped[['SampleID', 'Population']].copy()
    sample_map.columns = ['sample_id', 'pop']
    
    # Save as CSV
    sample_map.to_csv(output_path, index=False)
    
    print(f"\n✅ Created sample map: {output_path}")
    print(f"   Total samples: {len(sample_map)}")
    print(f"   Populations: {len(sample_map['pop'].unique())}")
    print(f"\n   Population breakdown:")
    pop_counts = sample_map['pop'].value_counts().sort_index()
    for pop, count in pop_counts.items():
        print(f"     {pop}: {count} samples")
    
    # Show first few rows
    print(f"\n   First 5 rows:")
    print(sample_map.head().to_string(index=False))
    
    return output_path


if __name__ == "__main__":
    create_sample_map()

