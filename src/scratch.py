import pandas as pd

# Load patient mapping CSV
mapping_df = pd.read_csv("data/GRN/Astrocytes/EC_Astrocytes_patient_mapping.csv")

# Count number of unique patients
num_patients = mapping_df['patient_name']
print(f"Number of patients: {num_patients}")