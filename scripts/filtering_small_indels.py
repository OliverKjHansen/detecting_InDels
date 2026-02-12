import pandas as pd
import argparse
import sys

def load_file(filepath):
    """Load TSV file with proper type handling"""
    try:
        # df = pd.read_csv(filepath, sep='\t')
        # return df
        df = pd.read_csv(filepath, sep='\t', # just a quick fix. remove later
                        usecols=range(29))
        print(f"Loaded {len(df.columns)} columns: {list(df.columns)}", file=sys.stderr)
        return df
    except FileNotFoundError:
        print(f"Error: File {filepath} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

def filter_indels(df, **filters):
    """
    Filter indels based on column conditions
    
    Parameters:
    df: DataFrame
    **filters: column_name=value or column_name=(operator, value)
              operators: 'eq', 'ne', 'gt', 'gte', 'lt', 'lte', 'in', 'notin'
    """
    filtered_df = df.copy()
    
    for column, condition in filters.items():
        if column not in filtered_df.columns:
            print(f"Warning: Column '{column}' not found in data", file=sys.stderr)
            continue
        
        if isinstance(condition, tuple):
            operator, value = condition
            if operator == 'eq':
                filtered_df = filtered_df[filtered_df[column] == value]
            elif operator == 'ne':
                filtered_df = filtered_df[filtered_df[column] != value]
            elif operator == 'gt':
                filtered_df = filtered_df[filtered_df[column] > value]
            elif operator == 'gte':
                filtered_df = filtered_df[filtered_df[column] >= value]
            elif operator == 'lt':
                filtered_df = filtered_df[filtered_df[column] < value]
            elif operator == 'lte':
                filtered_df = filtered_df[filtered_df[column] <= value]
            elif operator == 'in':
                filtered_df = filtered_df[filtered_df[column].isin(value)]
            elif operator == 'notin':
                filtered_df = filtered_df[~filtered_df[column].isin(value)]
        else:
            # Simple equality
            filtered_df = filtered_df[filtered_df[column] == condition]
    
    return filtered_df

def main():
    parser = argparse.ArgumentParser(
        description='Filter small indels based on haploid mapping results'
    )
    parser.add_argument('indel_file', help='Small indels TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output filtered file')
    parser.add_argument('--haploid-check', choices=['True', 'False'], 
                        help='Filter by haploid_check column')
    parser.add_argument('--softclip', choices=['True', 'False'],
                        help='Filter by softclip_haploid column')
    parser.add_argument('--min-distance', type=int, default=0,
                        help='Minimum distance to indel in haploid (default: 0)')
    parser.add_argument('--max-distance', type=int,
                        help='Maximum distance to indel in haploid')
    parser.add_argument('--indel-type', choices=['insertion', 'deletion'],
                        help='Filter by indel type')
    parser.add_argument('--tissue', 
                        help='Filter by tissue (comma-separated for multiple)')
    parser.add_argument('--unique-only', action='store_true',
                        help='Keep only unique indels')
    parser.add_argument('--min-quality', type=int,
                        help='Minimum kmerQ quality score')
    
    args = parser.parse_args()
    
    # Load file
    print(f"Loading {args.indel_file}...", file=sys.stderr)
    df = load_file(args.indel_file)
    print(f"Loaded {len(df)} rows", file=sys.stderr)
    
    # Build filters
    filters = {}
    
    if args.haploid_check:
        filters['haploid_check'] = args.haploid_check == 'True'
    
    if args.softclip:
        filters['softclip_haploid'] = args.softclip == 'True'
    
    if args.min_distance is not None:
        filters['distance_to_indel_in_hap'] = ('gte', args.min_distance)
    
    if args.max_distance is not None:
        filters['distance_to_indel_in_hap'] = ('lte', args.max_distance)
    
    if args.indel_type:
        filters['indel_type'] = args.indel_type
    
    if args.tissue:
        tissues = [t.strip() for t in args.tissue.split(',')]
        filters['tissue'] = ('in', tissues)
    
    if args.unique_only:
        filters['unique'] = 'unique'
    
    if args.min_quality:
        filters['kmerQ'] = ('gte', args.min_quality)
    
    # Apply filters
    print(f"Applying {len(filters)} filters...", file=sys.stderr)
    filtered_df = filter_indels(df, **filters)
    print(f"Filtered to {len(filtered_df)} rows", file=sys.stderr)
    
    # Save output
    filtered_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()