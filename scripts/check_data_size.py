#!/usr/bin/env python3
"""
Quick script to check your data size and suggest how many more SNPs you need.
"""

import pandas as pd
import sys
import os

# Add src to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(os.path.join(project_root, 'src'))


def check_data_size(csv_path):
    """Check the size and coverage of your allele frequency data."""
    print("="*60)
    print("DATA SIZE CHECKER")
    print("="*60)
    
    df = pd.read_csv(csv_path)
    
    # Basic stats
    unique_positions = df['pos'].nunique()
    total_rows = len(df)
    populations = df['pop'].unique()
    region_start = df['pos'].min()
    region_end = df['pos'].max()
    region_size = region_end - region_start
    
    print(f"\n📊 Current Dataset:")
    print(f"  Unique SNP positions: {unique_positions}")
    print(f"  Total rows: {total_rows}")
    print(f"  Populations: {', '.join(populations)}")
    print(f"  Region: chr{df['chrom'].iloc[0]} {region_start:,} - {region_end:,}")
    print(f"  Region size: {region_size:,} bp ({region_size/1000:.1f} kb)")
    
    # Recommendations
    print(f"\n💡 Recommendations:")
    if unique_positions < 20:
        print(f"  ⚠️  You have {unique_positions} SNPs - this is quite small for HMM analysis")
        print(f"  📈 Target: 50-100 SNPs for robust results")
        print(f"  ➕ You need approximately {50 - unique_positions} - {100 - unique_positions} more SNPs")
    elif unique_positions < 50:
        print(f"  ✓ You have {unique_positions} SNPs - this is okay but could be better")
        print(f"  📈 Target: 50-100 SNPs for robust results")
        print(f"  ➕ You need approximately {50 - unique_positions} - {100 - unique_positions} more SNPs")
    elif unique_positions < 100:
        print(f"  ✓ You have {unique_positions} SNPs - this is good!")
        print(f"  📈 Target: 100+ SNPs for publication-quality results")
        print(f"  ➕ Consider adding {100 - unique_positions} more SNPs if possible")
    else:
        print(f"  ✓✓ You have {unique_positions} SNPs - excellent size for HMM analysis!")
    
    # Suggest region expansion
    print(f"\n🗺️  Suggested Region Expansion:")
    current_center = (region_start + region_end) / 2
    current_size = region_end - region_start
    
    # Expand by 2x
    expanded_start = int(current_center - current_size)
    expanded_end = int(current_center + current_size)
    print(f"  Current: {region_start:,} - {region_end:,} ({current_size/1000:.1f} kb)")
    print(f"  Expand 2x: {expanded_start:,} - {expanded_end:,} ({(expanded_end-expanded_start)/1000:.1f} kb)")
    
    # Expand by 5x (recommended)
    expanded_start_5x = int(current_center - current_size * 2.5)
    expanded_end_5x = int(current_center + current_size * 2.5)
    print(f"  Expand 5x: {expanded_start_5x:,} - {expanded_end_5x:,} ({(expanded_end_5x-expanded_start_5x)/1000:.1f} kb) ⭐ Recommended")
    
    # Check data quality
    print(f"\n🔍 Data Quality Check:")
    missing_af = df[df['af'].isna()].shape[0]
    if missing_af > 0:
        print(f"  ⚠️  {missing_af} rows have missing allele frequencies")
    else:
        print(f"  ✓ No missing allele frequencies")
    
    # Check if all populations have data for all positions
    pos_counts = df.groupby('pos')['pop'].count()
    expected_per_pos = len(populations)
    incomplete = (pos_counts < expected_per_pos).sum()
    if incomplete > 0:
        print(f"  ⚠️  {incomplete} positions are missing data for some populations")
    else:
        print(f"  ✓ All positions have data for all populations")
    
    print("\n" + "="*60)
    print("Next Steps:")
    print("1. Process your data: python src/data_prep.py --af-table data/example_abo_af.csv --pop-a YRI --pop-b CEU --out results/abo_freqs")
    print("2. Expand region to get more SNPs (see GET_MORE_DATA.md)")
    print("3. Load and analyze: from data_prep import load_hmm_data")
    print("="*60)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        csv_path = sys.argv[1]
    else:
        csv_path = os.path.join(project_root, 'data', 'example_abo_af.csv')
    
    if not os.path.exists(csv_path):
        print(f"Error: File not found: {csv_path}")
        sys.exit(1)
    
    check_data_size(csv_path)

