"""
Find and print metadata for Astrocytes clusters using Astrocytes metadata file.
"""

import pandas as pd
import os
import pickle

# Path to Astrocytes metadata file
METADATA_FILE = os.path.join("data", "AtlasData", "2025-11-16_Astrocytes_metadata.xlsx")
CACHE_FILE = os.path.join("data", "AtlasData", "2025-11-16_Astrocytes_metadata_cache.pkl")

def load_metadata(use_cache=True):
    """Load the Astrocytes metadata file (uses cache if available for faster loading)."""
    if not os.path.exists(METADATA_FILE):
        raise FileNotFoundError(f"Metadata file not found: {METADATA_FILE}")
    
    # Check if cache exists and is newer than source file
    use_cache_file = False
    if use_cache and os.path.exists(CACHE_FILE):
        cache_time = os.path.getmtime(CACHE_FILE)
        source_time = os.path.getmtime(METADATA_FILE)
        if cache_time > source_time:
            use_cache_file = True
            print(f"Loading metadata from cache: {CACHE_FILE}")
    
    if use_cache_file:
        # Load from cache (much faster)
        metadata_df = pd.read_pickle(CACHE_FILE)
    else:
        # Load from Excel and save to cache
        print(f"Loading metadata from Excel (this may take a while)...")
        metadata_df = pd.read_excel(METADATA_FILE, sheet_name=0)
        print(f"Saving to cache: {CACHE_FILE}")
        os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
        pd.to_pickle(metadata_df, CACHE_FILE)
        print(f"Cache saved. Future runs will be faster.")
    
    return metadata_df

def load_cluster_data(cluster_id):
    """Load the cluster metadata from pickle file."""
    cluster_pickle = f"data/Astrocyte_cluster_{cluster_id}_cells.pkl"
    if not os.path.exists(cluster_pickle):
        raise FileNotFoundError(f"Cluster pickle file not found: {cluster_pickle}")
    
    with open(cluster_pickle, 'rb') as f:
        cluster_data = pickle.load(f)
    return cluster_data

def find_braak_stage(data, four_digit_ids):
    """Find the Braak stage distribution of the cells that are in the cluster."""
    braak_stages = {}
    
    if 'cell_annotation' not in data.columns or 'Braak' not in data.columns:
        print("Warning: Required columns not found")
        return braak_stages
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            braak_value = matching_rows['Braak'].iloc[0]
            if pd.notna(braak_value):
                braak_str = str(braak_value).strip()
                # Initialize stage if not seen before, then increment
                if braak_str not in braak_stages:
                    braak_stages[braak_str] = 0
                braak_stages[braak_str] += 1
    
    return braak_stages

def find_thal_distribution(data, four_digit_ids):
    """Find the Thal distribution of the cells that are in the cluster."""
    thaal_distribution = {}
    if 'cell_annotation' not in data.columns or 'Thal' not in data.columns:
        print("Warning: Required columns not found")
        return thaal_distribution
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            thaal_value = matching_rows['Thal'].iloc[0]
            if pd.notna(thaal_value):
                thaal_str = str(thaal_value).strip()
                if thaal_str not in thaal_distribution:
                    thaal_distribution[thaal_str] = 0
                thaal_distribution[thaal_str] += 1
    return thaal_distribution

def find_patient_distribution(data, four_digit_ids):
    """Find the patient distribution of the cells that are in the cluster."""
    patient_distribution = {}
    if 'cell_annotation' not in data.columns:
        print("Warning: Required columns not found")
        return patient_distribution
    
    # Try to find patient column (check common names)
    patient_column = 'Donor.ID.y'
    
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            patient_value = matching_rows[patient_column].iloc[0]
            if pd.notna(patient_value):
                patient_str = str(patient_value).strip()
                if patient_str not in patient_distribution:
                    patient_distribution[patient_str] = 0
                patient_distribution[patient_str] += 1
    return patient_distribution

def find_stage_distribution(data, four_digit_ids):
    """Find the Pathology Stage distribution of the cells that are in the cluster."""
    stage_distribution = {}
    if 'cell_annotation' not in data.columns:
        print("Warning: Required columns not found")
        return stage_distribution
    
    # Try to find stage column (check common names)
    stage_column = 'Path..Group.'
    
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            stage_value = matching_rows[stage_column].iloc[0]
            if pd.notna(stage_value):
                stage_str = str(stage_value).strip()
                if stage_str not in stage_distribution:
                    stage_distribution[stage_str] = 0
                stage_distribution[stage_str] += 1
    return stage_distribution

def find_total_tau_distribution(data, four_digit_ids):
    """Find the mean and SD of TOTAL.TAU..ng.mg. values for the given four-digit IDs."""
    import numpy as np
    
    tau_values = []
    if 'cell_annotation' not in data.columns or 'TOTAL.TAU..ng.mg.' not in data.columns:
        if 'TOTAL.TAU..ng.mg.' not in data.columns:
            print("Warning: 'TOTAL.TAU..ng.mg.' column not found")
        return None, None
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            tau_value = matching_rows['TOTAL.TAU..ng.mg.'].iloc[0]
            if pd.notna(tau_value):
                tau_values.append(float(tau_value))
    
    if len(tau_values) == 0:
        return None, None
    
    mean_val = np.mean(tau_values)
    std_val = np.std(tau_values, ddof=1) if len(tau_values) > 1 else 0.0
    
    return mean_val, std_val

def find_ptau_total_tau_distribution(data, four_digit_ids):
    """Find the mean and SD of Ptau.Total.Tau..A.U. values for the given four-digit IDs."""
    import numpy as np
    
    ptau_values = []
    if 'cell_annotation' not in data.columns or 'Ptau.Total.Tau..A.U..' not in data.columns:
        if 'Ptau.Total.Tau..A.U.' not in data.columns:
            print("Warning: 'Ptau.Total.Tau..A.U..' column not found")
        return None, None
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            ptau_value = matching_rows['Ptau.Total.Tau..A.U..'].iloc[0]
            if pd.notna(ptau_value):
                ptau_values.append(float(ptau_value))
    
    if len(ptau_values) == 0:
        return None, None
    
    mean_val = np.mean(ptau_values)
    std_val = np.std(ptau_values, ddof=1) if len(ptau_values) > 1 else 0.0
    
    return mean_val, std_val

def find_ad_distribution(data, four_digit_ids):
    """Find the AD distribution of the cells that are in the cluster."""
    ad_distribution = {}
    if 'cell_annotation' not in data.columns or 'AD' not in data.columns:
        print("Warning: Required columns not found")
        return ad_distribution
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            ad_value = matching_rows['AD'].iloc[0]
            if pd.notna(ad_value):
                ad_str = str(ad_value).strip()
                if ad_str not in ad_distribution:
                    ad_distribution[ad_str] = 0
                ad_distribution[ad_str] += 1
    return ad_distribution

def find_apoe_distribution(data, four_digit_ids):
    """Find the APOE (E column) distribution of the cells that are in the cluster."""
    apoe_distribution = {}
    if 'cell_annotation' not in data.columns or 'E' not in data.columns:
        print("Warning: Required columns not found")
        return apoe_distribution
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            apoe_value = matching_rows['E'].iloc[0]
            if pd.notna(apoe_value):
                apoe_str = str(apoe_value).strip()
                if apoe_str not in apoe_distribution:
                    apoe_distribution[apoe_str] = 0
                apoe_distribution[apoe_str] += 1
    return apoe_distribution

def find_age_distribution(data, four_digit_ids):
    """Find the mean and SD of age values for the given four-digit IDs."""
    import numpy as np
    
    age_values = []
    # Try to find age column (check common names)
    age_column = None
    for col in ['Age', 'age', 'Age.at.Death', 'Age_at_Death']:
        if col in data.columns:
            age_column = col
            break
    
    if age_column is None or 'cell_annotation' not in data.columns:
        if age_column is None:
            print("Warning: Age column not found")
        return None, None
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            age_value = matching_rows[age_column].iloc[0]
            if pd.notna(age_value):
                # Handle values like '90+' by removing the '+' and converting
                age_str = str(age_value).strip().replace('+', '')
                try:
                    age_values.append(float(age_str))
                except ValueError:
                    # Skip values that can't be converted to float
                    continue
    
    if len(age_values) == 0:
        return None, None
    
    mean_val = np.mean(age_values)
    std_val = np.std(age_values, ddof=1) if len(age_values) > 1 else 0.0
    
    return mean_val, std_val

def find_sex_distribution(data, four_digit_ids):
    """Find the sex distribution of the cells that are in the cluster."""
    sex_distribution = {}
    if 'cell_annotation' not in data.columns:
        print("Warning: Required columns not found")
        return sex_distribution
    
    # Try to find sex column (check common names)
    sex_column = None
    for col in ['Sex', 'sex', 'Gender', 'gender']:
        if col in data.columns:
            sex_column = col
            break
    
    if sex_column is None:
        print("Warning: Sex column not found")
        return sex_distribution
    
    for four_digit_id in four_digit_ids:
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            sex_value = matching_rows[sex_column].iloc[0]
            if pd.notna(sex_value):
                sex_str = str(sex_value).strip()
                if sex_str not in sex_distribution:
                    sex_distribution[sex_str] = 0
                sex_distribution[sex_str] += 1
    return sex_distribution

def find_pmi_distribution(data, four_digit_ids):
    """Find the mean and SD of PMI values for the given four-digit IDs."""
    import numpy as np
    
    pmi_values = []
    # Try to find PMI column (check common names)
    pmi_column = None
    for col in ['PMI', 'pmi', 'Post.Mortem.Interval', 'Post_Mortem_Interval']:
        if col in data.columns:
            pmi_column = col
            break
    
    if pmi_column is None or 'cell_annotation' not in data.columns:
        if pmi_column is None:
            print("Warning: PMI column not found")
        return None, None
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            pmi_value = matching_rows[pmi_column].iloc[0]
            if pd.notna(pmi_value):
                try:
                    pmi_values.append(float(pmi_value))
                except ValueError:
                    # Skip values that can't be converted to float
                    continue
    
    if len(pmi_values) == 0:
        return None, None
    
    mean_val = np.mean(pmi_values)
    std_val = np.std(pmi_values, ddof=1) if len(pmi_values) > 1 else 0.0
    
    return mean_val, std_val

def find_rin_distribution(data, four_digit_ids):
    """Find the mean and SD of RIN values for the given four-digit IDs."""
    import numpy as np
    
    rin_values = []
    # Try to find RIN column (check common names)
    rin_column = None
    for col in ['RIN', 'rin', 'RNA.Integrity.Number', 'RNA_Integrity_Number']:
        if col in data.columns:
            rin_column = col
            break
    
    if rin_column is None or 'cell_annotation' not in data.columns:
        if rin_column is None:
            print("Warning: RIN column not found")
        return None, None
    
    for four_digit_id in four_digit_ids:
        # Match cell_annotation that contains this 4-digit ID
        matching_rows = data[data['cell_annotation'].astype(str).str.contains(four_digit_id, na=False)]
        if len(matching_rows) > 0:
            rin_value = matching_rows[rin_column].iloc[0]
            if pd.notna(rin_value):
                try:
                    rin_values.append(float(rin_value))
                except ValueError:
                    # Skip values that can't be converted to float
                    continue
    
    if len(rin_values) == 0:
        return None, None
    
    mean_val = np.mean(rin_values)
    std_val = np.std(rin_values, ddof=1) if len(rin_values) > 1 else 0.0
    
    return mean_val, std_val

if __name__ == "__main__":
    import re
    
    metadata = load_metadata()
    for cluster_id in range(7):
        cluster_data = load_cluster_data(cluster_id)
        
        # Extract 4-digit IDs from cell names once (format: "6289-MW-0053_bin2" -> "0053")
        # Pattern: -XXXX_ where XXXX is 4 digits (matches "-0053_" in "6289-MW-0053_bin2")
        cell_names = cluster_data.get('cell_names', [])
        four_digit_ids = []
        for cell_name in cell_names:
            # Match pattern like -0053_ or -0053-
            match = re.search(r'-(\d{4})[_-]', str(cell_name))
            if match:
                four_digit_ids.append(match.group(1))
        
        mean_rin, std_rin = find_rin_distribution(metadata, four_digit_ids)
        if mean_rin is not None:
            print(f"Cluster {cluster_id}: RIN mean={mean_rin:.2f}, SD={std_rin:.2f}")
        else:
            print(f"Cluster {cluster_id}: No RIN values found")