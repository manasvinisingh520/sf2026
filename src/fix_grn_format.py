"""
Fix GRN file format to ensure each row is properly tab-separated.

The GRN files currently have TF on one line and Target/Score on the next line.
This script converts them to proper tab-separated format: TF\tTarget\tScore
"""

import os
import argparse
from pathlib import Path
from utils import get_grn_file_paths


def fix_grn_format(input_file, output_file):
    """
    Fix GRN file format to ensure tab-separated rows.
    Pattern: TF on one line, Target Score on next line (with leading space)
    Converts to: TF\tTarget\tScore
    
    Parameters:
    -----------
    input_file : str
        Path to input GRN file
    output_file : str
        Path to output GRN file
    """
    print(f"\nReading GRN file: {input_file}")
    file_size = os.path.getsize(input_file) / (1024 * 1024)
    print(f"File size: {file_size:.2f} MB")
    
    # Count TFs in original file
    tf_count = 0
    rows_written = 0
    
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8') as outfile:
        
        while True:
            # Read TF line
            tf_line = infile.readline()
            if not tf_line:  # End of file
                break
            
            tf_line = tf_line.rstrip('\n\r')
            if not tf_line:  # Empty line, skip
                continue
            
            # Check if already tab-separated
            if '\t' in tf_line:
                # Already in correct format, write as is
                outfile.write(tf_line + '\n')
                rows_written += 1
                continue
            
            # Read Target/Score line
            target_line = infile.readline()
            if not target_line:  # End of file, no Target/Score
                break
            
            target_line = target_line.rstrip('\n\r')
            if not target_line:  # Empty line, skip
                continue
            
            # Remove leading whitespace and split
            parts = target_line.strip().split(None, 1)
            if len(parts) == 2:
                target, score = parts
                # Write as tab-separated: TF\tTarget\tScore
                outfile.write(f"{tf_line}\t{target}\t{score}\n")
                rows_written += 1
                tf_count += 1
            else:
                print(f"  Warning: Target/Score line '{target_line}' at line {rows_written + 1} has no Score")
                quit()
    
    print(f"Processed {rows_written} rows")
    print(f"Successfully wrote fixed format to: {output_file}")
    
    # Validation check
    print(f"\nValidation:")
    print(f"  TFs found in original file: {tf_count}")
    print(f"  Rows written to new file: {rows_written}")
    if tf_count == rows_written:
        print(f"  ✓ Match! All TFs converted successfully.")
    else:
        print(f"  ⚠ Mismatch: {abs(tf_count - rows_written)} difference")


def main():
    """Main function to fix GRN file format."""
    parser = argparse.ArgumentParser(
        description='Fix GRN file format to ensure proper tab-separated rows',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-region',
        type=str,
        required=True,
        help='Region to process (e.g., V1, EC, ITG, PFC, V2)'
    )
    parser.add_argument(
        '-cell_type',
        type=str,
        default='Astrocytes',
        help='Cell type (e.g., Astrocytes, Microglia). Default: Astrocytes'
    )
    parser.add_argument(
        '-input_dir',
        type=str,
        default=None,
        help='Input directory (default: data)'
    )
    parser.add_argument(
        '-output_suffix',
        type=str,
        default='_fixed',
        help='Suffix to add to output filename (default: _fixed). Use empty string to overwrite original.'
    )
    
    args = parser.parse_args()
    
    # Determine data directory
    if args.input_dir:
        data_dir = args.input_dir
    else:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    
    # Get file paths using utility function
    grn_file, _, _, _, _, _ = get_grn_file_paths(
        args.region, cell_type=args.cell_type, data_dir=data_dir
    )
    
    # Create output filename
    if args.output_suffix:
        # Add suffix before extension
        base, ext = os.path.splitext(grn_file)
        output_file = f"{base}{args.output_suffix}{ext}"
    else:
        # Overwrite original
        output_file = grn_file
    
    if not os.path.exists(grn_file):
        print(f"Error: GRN file not found: {grn_file}")
        return
    
    try:
        fix_grn_format(grn_file, output_file)
        print("\n" + "="*60)
        print("Complete!")
        print("="*60)
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
