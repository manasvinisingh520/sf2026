import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import re


def get_top_k_genes(df: pd.DataFrame, k: int, sort_by: str = 'padj') -> set:
    df_clean = df[df[sort_by].notna()].copy()
    if len(df_clean) == 0:
        return set()
    ascending = (sort_by == 'padj')
    df_sorted = df_clean.sort_values(sort_by, ascending=ascending)
    return set(df_sorted.head(k).index)


def analyze_bin_overlap(results_dir: str = "results/dge4", top_k: int = 50):
    results_dir = Path(results_dir)
    if not results_dir.exists():
        raise ValueError(f"Directory {results_dir} does not exist")
    
    csv_files = sorted(list(results_dir.glob("*.csv")))
    if len(csv_files) == 0:
        raise ValueError(f"No CSV files found in {results_dir}")
    
    file_groups = defaultdict(list)
    for file_path in csv_files:
        filename = file_path.stem.replace('dge_results_', '')
        match = re.match(r'(.+?)_bins(\d+)(?:_seed(\d+))?$', filename)
        if match:
            comparison = match.group(1)
            bin_size = int(match.group(2))
            seed = match.group(3) if match.group(3) else "100"
            file_groups[(comparison, bin_size)].append((file_path, seed))
    
    results = {}
    for (comparison, bin_size), files in file_groups.items():
        if len(files) < 2:
            continue
        
        seed_top_genes = {}
        for file_path, seed in sorted(files, key=lambda x: x[1]):
            try:
                df = pd.read_csv(file_path, index_col=0)
                if 'padj' not in df.columns:
                    continue
                seed_top_genes[seed] = get_top_k_genes(df, top_k, sort_by='padj')
            except Exception:
                continue
        
        if len(seed_top_genes) < 2:
            continue
        
        common_genes = set.intersection(*seed_top_genes.values())
        overlap_pct = (len(common_genes) / top_k) * 100
        
        print(f"{comparison}, bins={bin_size}: {overlap_pct:.2f}%")
        
        results[(comparison, bin_size)] = {
            'n_seeds': len(seed_top_genes),
            'overlap_pct': overlap_pct
        }
    
    return results


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze gene overlap across seeds for different bin sizes')
    parser.add_argument('--results_dir', type=str, default='results/dge4',
                       help='Directory containing DGE result CSV files (default: results/dge4)')
    parser.add_argument('--top_k', type=int, default=50,
                       help='Number of top genes to consider (default: 50)')
    
    args = parser.parse_args()
    results = analyze_bin_overlap(args.results_dir, args.top_k)
    
    # Calculate summary statistics by bin size
    bin_summary = defaultdict(list)
    for (comparison, bin_size), data in results.items():
        bin_summary[bin_size].append(data['overlap_pct'])
    
    print()
    for bin_size in sorted(bin_summary.keys()):
        overlap_pcts = bin_summary[bin_size]
        mean_overlap = np.mean(overlap_pcts)
        std_overlap = np.std(overlap_pcts)
        print(f"bins={bin_size}: mean={mean_overlap:.2f}%, std={std_overlap:.2f}% (n={len(overlap_pcts)})")


if __name__ == "__main__":
    main()

