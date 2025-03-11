import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# Example Dataset (Replace with your actual dataset)
data =pd.read_csv(r"C:\Users\ashir\Music\updated_ddi_dataset_100.csv")

# Convert to DataFrame
df = pd.DataFrame(data)

# Function to convert SMILES to Morgan fingerprints
def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, n_bits))

# Feature Engineering: Combine fingerprints of two drugs
def combine_fingerprints(smiles1, smiles2):
    fp1 = smiles_to_fingerprint(smiles1)
    fp2 = smiles_to_fingerprint(smiles2)
    return np.concatenate((fp1, fp2))

# Generate features and labels
X = np.array([combine_fingerprints(d1, d2) for d1, d2 in zip(df["smiles1"], df["smiles2"])])
y = df["Side Effect Name"].apply(lambda x: pd.Series({effect: 1 for effect in x})).fillna(0).values

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a multi-label classifier
classifier = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, random_state=42))
classifier.fit(X_train, y_train)

# Evaluate the model
y_pred = classifier.predict(X_test)
print(classification_report(y_test, y_pred))

# Function to predict side effects for new drug pairs
def predict_side_effects(smiles1, smiles2):
    features = combine_fingerprints(smiles1, smiles2).reshape(1, -1)
    predictions = classifier.predict(features)
    side_effects = df["Side Effect Name"].explode().unique()
    return [effect for effect, pred in zip(side_effects, predictions[0]) if pred == 1]

