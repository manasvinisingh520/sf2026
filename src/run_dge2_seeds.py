"""
Run dge2.py for multiple seeds (5-49).

This script calls dge2.py repeatedly with different seed values.
"""

import subprocess
import sys
import argparse
from pathlib import Path


def run_dge2_for_seeds(group1, group2, region, seed_start=5, seed_end=49, 
                       bins=800):
    """
    Run dge2.py for seeds in the specified range.
    
    Parameters:
    -----------
    group1 : str
        First group name for comparison
    group2 : str
        Second group name for comparison
    region : str
        Region identifier (e.g., EC, ITG, PFC, V1, V2)
    seed_start : int
        Starting seed value (default: 5)
    seed_end : int
        Ending seed value (default: 49, inclusive)
    bins : int
        Target number of cells per bin (default: 800)
    """
    script_path = Path(__file__).parent / "dge2.py"
    
    if not script_path.exists():
        raise FileNotFoundError(f"dge2.py not found at {script_path}")
    
    seeds = list(range(seed_start, seed_end + 1))
    total_seeds = len(seeds)
    
    print(f"Running dge2.py for seeds {seed_start}-{seed_end} ({total_seeds} seeds)")
    print(f"Parameters: group1={group1}, group2={group2}, region={region}, bins={bins}")
    print()
    
    failed_seeds = []
    
    for i, seed in enumerate(seeds, 1):
        print(f"[{i}/{total_seeds}] Running seed {seed}...")
        
        cmd = [
            sys.executable,
            str(script_path),
            "--group1", str(group1),
            "--group2", str(group2),
            "--region", str(region),
            "--seed", str(seed),
            "--bins", str(bins)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            print(f"  ✓ Seed {seed} completed successfully")
            if result.stdout:
                # Print last few lines of output if available
                lines = result.stdout.strip().split('\n')
                if lines:
                    print(f"    {lines[-1]}")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Seed {seed} failed with error:")
            print(f"    {e.stderr}")
            failed_seeds.append(seed)
        except Exception as e:
            print(f"  ✗ Seed {seed} failed with exception: {e}")
            failed_seeds.append(seed)
        
        print()
    
    # Summary
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total seeds: {total_seeds}")
    print(f"Successful: {total_seeds - len(failed_seeds)}")
    print(f"Failed: {len(failed_seeds)}")
    if failed_seeds:
        print(f"Failed seeds: {failed_seeds}")


def main():
    parser = argparse.ArgumentParser(
        description='Run dge2.py for multiple seeds (5-49)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        '--group1',
        type=str,
        required=True,
        help='First group name for comparison (e.g., "1")'
    )
    parser.add_argument(
        '--group2',
        type=str,
        required=True,
        help='Second group name for comparison (e.g., "2")'
    )
    parser.add_argument(
        '--region',
        type=str,
        required=True,
        help='Region identifier (e.g., EC, ITG, PFC, V1, V2)'
    )
    parser.add_argument(
        '--seed_start',
        type=int,
        default=5,
        help='Starting seed value (default: 5)'
    )
    parser.add_argument(
        '--seed_end',
        type=int,
        default=49,
        help='Ending seed value (default: 49, inclusive)'
    )
    parser.add_argument(
        '--bins',
        type=int,
        default=800,
        help='Target number of cells per bin for pseudobulk aggregation (default: 800)'
    )
    
    args = parser.parse_args()
    
    run_dge2_for_seeds(
        group1=args.group1,
        group2=args.group2,
        region=args.region,
        seed_start=args.seed_start,
        seed_end=args.seed_end,
        bins=args.bins
    )


if __name__ == "__main__":
    main()

