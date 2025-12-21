
import pandas as pd
import sys
import os


def filter_nonzero_deltaaf(input_path, output_path=None):

    print(f"Loading data from: {input_path}")
    df = pd.read_csv(input_path)
    
    print(f"Original data: {len(df)} SNPs")
    print(f"  Zero DeltaAF: {(df['delta_af'] == 0).sum()} SNPs")
    print(f"  Non-zero DeltaAF: {(df['delta_af'] > 0).sum()} SNPs")
    
    filtered_df = df[df['delta_af'] > 0].copy()
    
    print(f"\nFiltered data: {len(filtered_df)} SNPs (non-zero only)")
    print(f"  DeltaAF range: {filtered_df['delta_af'].min():.6f} - {filtered_df['delta_af'].max():.6f}")
    print(f"  Mean DeltaAF: {filtered_df['delta_af'].mean():.6f}")
    
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_nonzero{ext}"
    
    filtered_df.to_csv(output_path, index=False)
    print(f"\n✅ Saved filtered data to: {output_path}")
    
    return filtered_df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Filter delta_af.csv to non-zero values")
    parser.add_argument('--input', type=str, 
                       default='results/abo_expanded/delta_af.csv',
                       help='Input delta_af.csv file')
    parser.add_argument('--output', type=str, default=None,
                       help='Output file path (default: adds _nonzero suffix)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ Error: {args.input} not found!")
        sys.exit(1)
    
    filter_nonzero_deltaaf(args.input, args.output)

