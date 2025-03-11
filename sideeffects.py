import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import pickle

class SideEffectPredictor:
    def __init__(self):
        self.mlb = MultiLabelBinarizer()
        self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        
    def generate_morgan_fingerprints(self, smiles):
        """Convert SMILES to Morgan fingerprints"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))
    
    def prepare_data(self, df):
        """Prepare the data for training"""
        # Group by SMILES and aggregate side effects
        grouped = df.groupby('SMILES')['Side_Effect_Name'].agg(list).reset_index()
        
        # Generate fingerprints for each SMILES
        X = np.array([self.generate_morgan_fingerprints(s) for s in grouped['SMILES']])
        
        # Transform side effects to multi-label format
        y = self.mlb.fit_transform(grouped['Side_Effect_Name'])
        
        return X, y
    
    def train(self, df):
        """Train the model"""
        X, y = self.prepare_data(df)
        
        # Split the data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        # Train the model
        self.model.fit(X_train, y_train)
        
        # Calculate and print accuracy
        train_score = self.model.score(X_train, y_train)
        test_score = self.model.score(X_test, y_test)
        print(f"Training accuracy: {train_score:.3f}")
        print(f"Testing accuracy: {test_score:.3f}")
        
    def predict(self, smiles):
        """Predict side effects for a given SMILES string"""
        # Generate fingerprints
        fp = self.generate_morgan_fingerprints(smiles)
        if fp is None:
            return "Invalid SMILES string"
            
        # Make prediction
        pred = self.model.predict([fp])[0]
        
        # Convert prediction back to side effect names
        side_effects = self.mlb.inverse_transform(pred.reshape(1, -1))[0]
        
        return list(side_effects)
    
    def save_model(self, filepath):
        """Save the trained model"""
        with open(filepath, 'wb') as f:
            pickle.dump({'model': self.model, 'mlb': self.mlb}, f)
    
    @classmethod
    def load_model(cls, filepath):
        """Load a trained model"""
        predictor = cls()
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
            predictor.model = data['model']
            predictor.mlb = data['mlb']
        return predictor

# Example usage:
def main():
    # Read the CSV file
    df = pd.read_csv('sider_smiles.csv')
    
    # Create and train the model
    predictor = SideEffectPredictor()
    predictor.train(df)
    
    # Save the model
    predictor.save_model('side_effect_model.pkl')
    
    # Example prediction
    test_smiles = input("Enter SMILE Structure: ")
    predicted_effects = predictor.predict(test_smiles)
    print(f"\nPredicted side effects for test compound:")
    for effect in predicted_effects:
        print(f"- {effect}")

if __name__ == "__main__":
    main()