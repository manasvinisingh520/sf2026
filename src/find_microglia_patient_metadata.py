"""
Find and print metadata for a list of patients using the combined Microglia metadata file.
"""

import pandas as pd
import os
import pickle

# Path to combined Microglia metadata file
METADATA_FILE = os.path.join("data", "AtlasData", "2025-12-27_Microglia_Metadata.csv")

def load_metadata():
    """Load the combined Microglia metadata file."""
    if not os.path.exists(METADATA_FILE):
        raise FileNotFoundError(f"Metadata file not found: {METADATA_FILE}")
    
    metadata_df = pd.read_csv(METADATA_FILE)
    return metadata_df

def load_cluster_metadata(cluster_id):
    """Load the cluster patients metadata from pickle file."""
    cluster_pickle = f"data/microglia_cluster_{cluster_id}_cells.pkl"
    with open(cluster_pickle, 'rb') as f:
        cluster_data = pickle.load(f)
    return cluster_data

def find_gender_distribution(data, cluster_data):
    male_count = 0
    female_count = 0
    unknown_count = 0
    
    # Get cell names from pickle (format: "2259_bin16" or "patient_number_binXX")
    patient_ids = cluster_data['cell_names']

    # Match patient numbers to Donor ID column
    # Convert Donor ID to string for matching
    data['Donor ID'] = data['Donor ID'].astype(str)
    
    # Count gender for each unique patient number
    for patient_number in patient_ids:
        if patient_number in data['Donor ID'].values:
            sex_value = data[data['Donor ID'] == patient_number]['Sex'].iloc[0]
            if sex_value == 'F':
                female_count += 1
            elif sex_value == 'M':
                male_count += 1
            else:
                unknown_count += 1
        else:
            unknown_count += 1 
    
    return male_count, female_count, unknown_count

def find_apoe_distribution(data, cluster_data):    
    patient_ids = cluster_data['cell_names']
    patient_ids = [patient_id.split('_')[0] for patient_id in patient_ids]
    data['Donor ID'] = data['Donor ID'].astype(str)
    apoe_23 = 0
    apoe_24 = 0
    apoe_34 = 0
    apoe_22 = 0
    apoe_33 = 0
    apoe_44 = 0
    for patient_number in patient_ids: 
        if patient_number in data['Donor ID'].values:
            apoe_value = data[data['Donor ID'] == patient_number]['APOE'].iloc[0]
            if apoe_value == '2/3':
                apoe_23 += 1
            elif apoe_value == '2/4':
                apoe_24 += 1
            elif apoe_value == '3/4':
                apoe_34 += 1
            elif apoe_value == '2/2':
                apoe_22 += 1
            elif apoe_value == '3/3':
                apoe_33 += 1
            elif apoe_value == '4/4':
                apoe_44 += 1
    
    return apoe_23, apoe_24, apoe_34, apoe_22, apoe_33, apoe_44

def find_stage_distribution(data, cluster_data):
    patient_ids = cluster_data['cell_names']
    patient_ids = [patient_id.split('_')[0] for patient_id in patient_ids]
    data['Donor ID'] = data['Donor ID'].astype(str)
    stage_count = {}
    stage_count['1'] = 0
    stage_count['2'] = 0
    stage_count['3'] = 0
    stage_count['4'] = 0
    for patient_number in patient_ids:
        if patient_number in data['Donor ID'].values:
            stage_value = str(data[data['Donor ID'] == patient_number][f'Pathology Stage'].iloc[0])
            stage_count[stage_value] += 1
    return stage_count

def find_braak_distribution(data, cluster_data):
    patient_ids = cluster_data['cell_names']
    patient_ids = [patient_id.split('_')[0] for patient_id in patient_ids]
    data['Donor ID'] = data['Donor ID'].astype(str)
    braak_count = {}
    braak_count['0'] = 0
    braak_count['I'] = 0
    braak_count['II'] = 0
    braak_count['III'] = 0
    braak_count['IV'] = 0
    braak_count['V'] = 0
    braak_count['VI'] = 0
    braak_count['unk'] = 0
    for patient_number in patient_ids:
        if patient_number in data['Donor ID'].values:
            braak_value = str(data[data['Donor ID'] == patient_number][f'Braak'].iloc[0])
            braak_count[braak_value] += 1
    return braak_count

def find_subclusters_distribution(data, cluster_data):
    patient_ids = cluster_data['cell_names']
    patient_ids = [patient_id.split('_')[0] for patient_id in patient_ids]
    data['Donor ID'] = data['Donor ID'].astype(str)
    subclusters_count = {}
    for patient_number in patient_ids:
        if patient_number in data['Donor ID'].values:
            subclusters_value = str(data[data['Donor ID'] == patient_number][f'Microglia Subclusters'].iloc[0])
            if subclusters_value not in subclusters_count:
                subclusters_count[subclusters_value] = 1
            else:
                subclusters_count[subclusters_value] += 1
    return subclusters_count

if __name__ == "__main__":
    metadata_df = load_metadata()
    for cluster_id in range(7):
        cluster_data = load_cluster_metadata(cluster_id)
        cell_names = cluster_data.get("cell_names", [])
        n_cells = len(cell_names)

        print(f"\n{'='*60}")
        print(f"Cluster {cluster_id} (n_cells={n_cells})")
        print(f"{'='*60}")

        male, female, unknown = find_gender_distribution(metadata_df, cluster_data)
        print(f"  Gender: M={male}, F={female}, unknown={unknown}")

        apoe_23, apoe_24, apoe_34, apoe_22, apoe_33, apoe_44 = find_apoe_distribution(metadata_df, cluster_data)
        apoe_dist = {"2/3": apoe_23, "2/4": apoe_24, "3/4": apoe_34, "2/2": apoe_22, "3/3": apoe_33, "4/4": apoe_44}
        print(f"  APOE: {apoe_dist}")

        stage = find_stage_distribution(metadata_df, cluster_data)
        print(f"  Stage (Pathology): {stage}")

        braak = find_braak_distribution(metadata_df, cluster_data)
        print(f"  Braak: {braak}")

        subclusters = find_subclusters_distribution(metadata_df, cluster_data)
        print(f"  Microglia subclusters: {subclusters}")