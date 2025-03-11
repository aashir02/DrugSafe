import pandas as pd
import pubchempy as pcp
from tqdm import tqdm  # Import tqdm for progress bar

# Function to get SMILES structure from STITCH ID
def get_smiles(stitch_id):
    try:
        if pd.isna(stitch_id) or not isinstance(stitch_id, str):
            return None  # Handle missing/invalid IDs
        
        cid = stitch_id.replace('CID', '')  # Remove 'CID' prefix
        compound = pcp.Compound.from_cid(cid)
        return compound.canonical_smiles
    except Exception as e:
        return None  # Return None in case of an error

# Load only the first 100 rows
df = pd.read_csv(r"C:\Users\ashir\Downloads\ChChSe-Decagon_polypharmacy.csv\ChChSe-Decagon_polypharmacy.csv", nrows=100)

# Apply function with progress bar
tqdm.pandas()  # Enable tqdm for pandas

print("Fetching SMILES for STITCH 1...")
df['smiles1'] = df['# STITCH 1'].progress_apply(get_smiles)

print("Fetching SMILES for STITCH 2...")
df['smiles2'] = df['STITCH 2'].progress_apply(get_smiles)

# Save the updated dataset
df.to_csv("updated_ddi_dataset_100.csv", index=False)

print("Dataset updated successfully with SMILES structures for 100 rows!")
