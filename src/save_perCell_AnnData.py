"""
Script to create, filter, and save AnnData object for a given region.
Supports Astrocytes (MTX + Excel metadata) and Microglia (MTX + per-region CSV metadata).
"""

import argparse
import os
import pandas as pd
from pathlib import Path
from utils import (
    read_mtx_file,
    create_anndata_object,
    read_excel_columns,
    get_region_file_paths,
)
from config import REGIONS, DATA_DIR, REGION_TO_TAB, METADATA_DATE_PREFIX

MICROGLIA_DATE_PREFIX = "2025-12-27"
MICROGLIA_REGIONS = ["CrossRegion", "ITG", "PFC"]


def _load_metadata_microglia(region, cell_names):
    """Load Microglia per-region CSV and align to MTX cell order."""
    metadata_path = os.path.join(DATA_DIR, f"{MICROGLIA_DATE_PREFIX}_Microglia_{region}_metadata.csv")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Microglia metadata not found: {metadata_path}")
    metadata = pd.read_csv(metadata_path)
    # Standardize column names
    rename = {}
    if "Donor ID" in metadata.columns:
        rename["Donor ID"] = "Donor_ID"
    if "Pathology Stage" in metadata.columns:
        rename["Pathology Stage"] = "Pathology_Stage"
    if rename:
        metadata = metadata.rename(columns=rename)
    # SampleName for compatibility with create_data / train_model
    if "SampleName" not in metadata.columns and "Donor_ID" in metadata.columns:
        metadata["SampleName"] = metadata["Donor_ID"].astype(str)
    # Thal / Ptau: map to names expected by create_data (Thal, Ptau.Total.Tau..A.U..)
    # Microglia may use CERAD (category) for Thal and Braak (category) for tau
    if "Thal" not in metadata.columns:
        for c in ["CERAD", "thal", "THAL", "Thal_phase"]:
            if c in metadata.columns:
                metadata["Thal"] = metadata[c].astype(str).str.strip()
                break
        if "Thal" not in metadata.columns:
            raise ValueError(f"Microglia metadata must have 'Thal' or 'CERAD'. Found: {list(metadata.columns)}")
    if "Ptau.Total.Tau..A.U.." not in metadata.columns:
        for c in ["Braak", "ptau", "Ptau", "Ptau.Total.Tau..A.U..", "Ptau_Total_Tau_A_U"]:
            if c not in metadata.columns:
                continue
            if c == "Braak":
                # Braak categories I, II, V, VI -> 2, 3, 5, 6 for regression
                braak_map = {"I": 2, "II": 3, "V": 5, "VI": 6}
                s = metadata["Braak"].astype(str).str.strip().str.upper()
                metadata["Ptau.Total.Tau..A.U.."] = s.replace(braak_map)
                metadata["Ptau.Total.Tau..A.U.."] = pd.to_numeric(metadata["Ptau.Total.Tau..A.U.."], errors="coerce")
                break
            else:
                metadata["Ptau.Total.Tau..A.U.."] = pd.to_numeric(metadata[c], errors="coerce")
                break
        if "Ptau.Total.Tau..A.U.." not in metadata.columns:
            raise ValueError(f"Microglia metadata must have 'Ptau.Total.Tau..A.U..' or 'Braak'. Found: {list(metadata.columns)}")
    # Align to cell_names: by cell_annotation column or by order
    cell_col = None
    for col in ["cell_annotation", "Cell", "Unnamed: 0", "cell"]:
        if col in metadata.columns and metadata[col].isin(cell_names).all():
            cell_col = col
            break
    if cell_col is not None:
        metadata = metadata.set_index(cell_col).reindex(cell_names)
    elif len(metadata) == len(cell_names):
        metadata = metadata.copy()
        metadata.index = cell_names
    else:
        raise ValueError(
            f"Microglia metadata rows ({len(metadata)}) must match cell count ({len(cell_names)}). "
            "Ensure CSV has a column matching MTX cell order or cell IDs."
        )
    # Filter: drop rows with missing Thal. For Astrocytes we also drop Thal == "unk".
    # Microglia CERAD has "none", "sparse", "frequent", "moderate" (all kept).
    valid = metadata["Thal"].notna() & (metadata["Thal"].astype(str).str.strip().str.lower() != "unk")
    metadata = metadata[valid].copy()
    # Ensure object/category columns are string so h5ad write doesn't fail (e.g. Braak mixed types)
    for col in metadata.columns:
        if metadata[col].dtype == object or isinstance(metadata[col].dtype, pd.CategoricalDtype):
            metadata[col] = metadata[col].astype(str)
    return metadata


def create_and_save_anndata(region, output_dir="data/model_data", date_prefix="2025-11-16", cell_type="Astrocytes"):

    print(f"\n{'='*60}")
    print(f"Processing region: {region} ({cell_type})")
    print(f"{'='*60}")

    is_microglia = cell_type == "Microglia"
    if is_microglia:
        date_prefix = MICROGLIA_DATE_PREFIX
        base_prefix = f"{date_prefix}_Microglia_{region}"
    else:
        base_prefix = f"{date_prefix}_Astrocytes_{region}"

    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=DATA_DIR,
        base_prefix=base_prefix
    )

    if is_microglia:
        # Load MTX first to get cell_names for metadata alignment
        matrix, gene_names, cell_names = read_mtx_file(
            mtx_path=str(mtx_path),
            row_annotation_path=str(row_annotation_path),
            col_annotation_path=str(col_annotation_path),
            transpose=False
        )
        print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
        filtered_metadata = _load_metadata_microglia(region, cell_names)
        # Filter matrix to valid cells (same order as filtered_metadata)
        valid_cell_names = filtered_metadata.index.tolist()
        cell_indices = [i for i, c in enumerate(cell_names) if c in valid_cell_names]
        cell_names = [cell_names[i] for i in cell_indices]
        matrix = matrix[:, cell_indices]
        print(f"Loaded metadata: {filtered_metadata.shape}, filtered to Thal != 'unk'")
    else:
        metadata_path = os.path.join(DATA_DIR, f"{date_prefix}_Astrocytes_metadata.xlsx")
        tab_index = REGION_TO_TAB.get(region)
        metadata = read_excel_columns(
            metadata_path,
            columns=['cell_annotation', 'Path..Group.', 'SampleName', 'percent.mito', 'RIN', 'Total.Genes.Detected', 'Thal', 'Ptau.Total.Tau..A.U..'],
            sheet_name=tab_index
        )
        print(f"Loaded metadata: {metadata.shape}")

        matrix, gene_names, cell_names = read_mtx_file(
            mtx_path=str(mtx_path),
            row_annotation_path=str(row_annotation_path),
            col_annotation_path=str(col_annotation_path),
            transpose=False
        )
        print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")

        valid_cell_mask = (metadata["Thal"] != "unk").fillna(False)
        valid_cell_names = metadata.loc[valid_cell_mask, "cell_annotation"].tolist()
        cell_indices = [i for i, c in enumerate(cell_names) if c in valid_cell_names]
        cell_names = [cell_names[i] for i in cell_indices]
        matrix = matrix[:, cell_indices]
        filtered_metadata = metadata[metadata["cell_annotation"].isin(cell_names)].set_index("cell_annotation").reindex(cell_names)

    # Keep gene names as from MTX (symbols). create_data converts to Ensembl when loading embeddings.
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,
        obs=filtered_metadata,
    )
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")

    print(f"AnnData object created: {adata.shape} (cells × genes)")

    output_filename = f"{region}_Microglia_AnnData_perCell.h5ad" if is_microglia else f"{region}_AnnData_perCell.h5ad"
    output_path = os.path.join(output_dir, output_filename)
    adata.write_h5ad(output_path)
    print(f"Successfully saved AnnData object to {output_path}")
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Create and save per-cell AnnData for a region (Astrocytes or Microglia)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-region", type=str, required=True, help="Region name (e.g. EC, ITG, CrossRegion for Microglia)")
    parser.add_argument("-cell_type", type=str, default="Astrocytes", choices=["Astrocytes", "Microglia"],
                        help="Cell type (default: Astrocytes)")
    parser.add_argument("-output_dir", type=str, default="data/model_data", help="Output directory for .h5ad file")
    args = parser.parse_args()

    region = args.region
    cell_type = args.cell_type
    if cell_type == "Astrocytes" and region not in REGIONS:
        parser.error(f"Astrocytes region must be one of {REGIONS}")
    if cell_type == "Microglia" and region not in MICROGLIA_REGIONS:
        parser.error(f"Microglia region must be one of {MICROGLIA_REGIONS}")

    try:
        create_and_save_anndata(
            region=region,
            output_dir=args.output_dir,
            cell_type=cell_type,
        )
    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
