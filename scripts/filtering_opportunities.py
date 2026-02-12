import pandas as pd
import argparse
import sys

def load_file(filepath):
    """Load TSV file with proper type handling"""
    try:
        df = pd.read_csv(filepath, sep='\t', low_memory=False,
                        on_bad_lines='skip')
        print(f"Loaded {len(df.columns)} columns: {list(df.columns)}", file=sys.stderr)
        return df
    except FileNotFoundError:
        print(f"Error: File {filepath} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

def filter_opportunities(df, **filters):
    """
    Filter opportunity data based on column conditions
    
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
        description='Filter opportunity file based on mapQ and kmerQ thresholds'
    )
    parser.add_argument('opportunity_file', help='Opportunity TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output filtered file')
    parser.add_argument('--min-mapq', type=int, default=0,
                        help='Minimum mapQ score (default: 0)')
    parser.add_argument('--max-mapq', type=int,
                        help='Maximum mapQ score')
    parser.add_argument('--min-kmerq', type=int, default=0,
                        help='Minimum kmerQ score (default: 0)')
    parser.add_argument('--max-kmerq', type=int,
                        help='Maximum kmerQ score')
    parser.add_argument('--tissue', 
                        help='Filter by tissue (comma-separated for multiple)')
    parser.add_argument('--context',
                        help='Filter by context (comma-separated for multiple)')
    parser.add_argument('--min-count', type=int,
                        help='Minimum count threshold')
    
    args = parser.parse_args()
    
    # Load file
    print(f"Loading {args.opportunity_file}...", file=sys.stderr)
    df = load_file(args.opportunity_file)
    print(f"Loaded {len(df)} rows", file=sys.stderr)
    
    # Build filters
    filters = {}
    
    if args.min_mapq is not None:
        filters['mapQ'] = ('gte', args.min_mapq)
    
    if args.max_mapq is not None:
        filters['mapQ'] = ('lte', args.max_mapq)
    
    if args.min_kmerq is not None:
        filters['kmerQ'] = ('gte', args.min_kmerq)
    
    if args.max_kmerq is not None:
        filters['kmerQ'] = ('lte', args.max_kmerq)
    
    if args.tissue:
        tissues = [t.strip() for t in args.tissue.split(',')]
        filters['tissue'] = ('in', tissues)
    
    if args.context:
        contexts = [c.strip() for c in args.context.split(',')]
        filters['context'] = ('in', contexts)
    
    if args.min_count:
        filters['count'] = ('gte', args.min_count)
    
    # Apply filters
    print(f"Applying {len(filters)} filters...", file=sys.stderr)
    filtered_df = filter_opportunities(df, **filters)
    print(f"Filtered to {len(filtered_df)} rows", file=sys.stderr)
    
    # Save output
    filtered_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()