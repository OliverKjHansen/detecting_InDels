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

def apply_filters_with_report(df, filter_steps):
    """
    Apply filters sequentially and report how many rows each step removes.
    
    Parameters:
        df: DataFrame
        filter_steps: list of (filter_name, column, condition) tuples
    
    Returns:
        filtered_df: DataFrame after all filters
        report: list of (filter_name, rows_before, rows_after, rows_removed) tuples
    """
    report = []
    current_df = df.copy()
    
    for filter_name, column, condition in filter_steps:
        rows_before = len(current_df)
        if column not in current_df.columns:
            print(f"Warning: Column '{column}' not found, skipping filter '{filter_name}'", file=sys.stderr)
            report.append((filter_name, rows_before, rows_before, 0))
            continue
        
        current_df = filter_indels(current_df, **{column: condition})
        rows_after = len(current_df)
        rows_removed = rows_before - rows_after
        report.append((filter_name, rows_before, rows_after, rows_removed))
        print(f"  {filter_name}: {rows_before} -> {rows_after} ({rows_removed} removed)", file=sys.stderr)
    
    return current_df, report

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
    parser.add_argument('--min-distance', type=int,
                        help='Minimum distance to indel in haploid')
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
    parser.add_argument('--min-read-end-distance', type=int,
                        help='Minimum distance from indel to read end (ToreadEnd column)')
    parser.add_argument('--report', 
                        help='Output filter report TSV file showing rows removed per step')
    
    args = parser.parse_args()
    
    # Load file
    print(f"Loading {args.indel_file}...", file=sys.stderr)
    df = load_file(args.indel_file)
    print(f"Loaded {len(df)} rows", file=sys.stderr)
    
    # Build filter steps as ordered list for reporting
    filter_steps = []
    
    if args.haploid_check:
        filter_steps.append(('haploid_check', 'haploid_check', args.haploid_check == 'True'))
    
    if args.softclip:
        filter_steps.append(('softclip_filter', 'softclip_haploid', args.softclip == 'True'))
    
    if args.min_distance is not None:
        filter_steps.append(('min_distance', 'distance_to_indel_in_hap', ('gte', args.min_distance)))
    
    if args.max_distance is not None:
        filter_steps.append(('max_distance', 'distance_to_indel_in_hap', ('lte', args.max_distance)))
    
    if args.indel_type:
        filter_steps.append(('indel_type', 'indel_type', args.indel_type))
    
    if args.tissue:
        tissues = [t.strip() for t in args.tissue.split(',')]
        filter_steps.append(('tissue', 'tissue', ('in', tissues)))
    
    if args.unique_only:
        filter_steps.append(('unique_only', 'unique', 'unique'))
    
    if args.min_quality:
        filter_steps.append(('min_kmerQ', 'kmerQ', ('gte', args.min_quality)))
    
    if args.min_read_end_distance:
        filter_steps.append(('min_read_end_distance', 'ToreadEnd', ('gte', args.min_read_end_distance)))
    
    # Apply filters with reporting
    print(f"Applying {len(filter_steps)} filters...", file=sys.stderr)
    filtered_df, report = apply_filters_with_report(df, filter_steps)
    print(f"Filtered to {len(filtered_df)} rows", file=sys.stderr)
    
    # Save output
    filtered_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved to {args.output}", file=sys.stderr)
    
    # Save filter report if requested
    if args.report:
        report_df = pd.DataFrame(report, columns=['filter', 'rows_before', 'rows_after', 'rows_removed'])
        report_df.to_csv(args.report, sep='\t', index=False)
        print(f"Filter report saved to {args.report}", file=sys.stderr)

if __name__ == "__main__":
    main()