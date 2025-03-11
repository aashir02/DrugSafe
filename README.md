
# Drug Side Effect & Interaction Predictor

This project provides a machine learning-based tool to predict side effects for single drugs and drug-drug interactions using SMILES (Simplified Molecular Input Line Entry System) structures. The tool is built using Python, Streamlit, and Scikit-learn, and it leverages Morgan fingerprints for molecular representation.

## Features

- **Single Drug Side Effect Prediction**: Predict side effects for a single drug using its SMILES structure.
- **Drug-Drug Interaction Prediction**: Predict side effects for combinations of two drugs using their SMILES structures.
- **User-Friendly Interface**: A Streamlit-based web application for easy interaction.
- **Machine Learning Models**: Utilizes Random Forest classifiers for predictions.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

### Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/drug-side-effect-predictor.git
   cd drug-side-effect-predictor
   ```

2. **Create a virtual environment** (optional but recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Download the pre-trained models**:
   - Place the pre-trained models (`side_effect_model.pkl` and any other required files) in the appropriate directory.
   - If you don't have the models, you can train them using the provided scripts (see [Training the Models](#training-the-models)).

## Usage

### Running the Streamlit App

1. **Start the Streamlit app**:
   ```bash
   streamlit run app.py
   ```

2. **Access the app**:
   - Open your web browser and go to `http://localhost:8501`.
   - Use the sidebar to navigate between **Single Drug Side Effects** and **Drug-Drug Interaction** prediction.

3. **Single Drug Side Effects**:
   - Enter a SMILES structure in the input box.
   - Click **Predict Side Effects** to see the predicted side effects.

4. **Drug-Drug Interaction**:
   - Enter SMILES structures for two drugs in the input boxes.
   - Click **Predict Interactions** to see the predicted side effects for the drug combination.

### Example SMILES Structures

Here are some example SMILES structures you can use for testing:

- **Aspirin**: `CC(=O)OC1=CC=CC=C1C(=O)O`
- **Paracetamol**: `CC(=O)NC1=CC=C(C=C1)O`
- **Ibuprofen**: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O`

### Training the Models

If you want to train the models from scratch, follow these steps:

1. **Prepare the dataset**:
   - Ensure you have the required datasets (`sider_smiles.csv` for single drug side effects and `updated_ddi_dataset_100.csv` for drug-drug interactions).
   - If you don't have the datasets, you can generate them using the `create_Dataset.py` script.

2. **Train the single drug side effect model**:
   - Run the `sideeffects.py` script to train the model:
     ```bash
     python sideeffects.py
     ```
   - This will generate the `side_effect_model.pkl` file.

3. **Train the drug-drug interaction model**:
   - Run the `Ddi_model.py` script to train the model:
     ```bash
     python Ddi_model.py
     ```
   - This will train the multi-label classifier for drug-drug interactions.

## Project Structure

```
drug-side-effect-predictor/
├── app.py                  # Streamlit application
├── sideeffects.py          # Single drug side effect predictor
├── Ddi_model.py            # Drug-drug interaction predictor
├── create_Dataset.py       # Script to generate SMILES dataset
├── requirements.txt        # Python dependencies
├── README.md               # Project documentation
├── side_effect_model.pkl   # Pre-trained single drug model
└── updated_ddi_dataset_100.csv  # Example dataset for drug-drug interactions
```

## Dependencies

The project relies on the following Python libraries:

- `streamlit`
- `pandas`
- `numpy`
- `rdkit`
- `scikit-learn`
- `pubchempy`
- `tqdm`

You can install all dependencies using:
```bash
pip install -r requirements.txt
```

## Contributing

Contributions are welcome! If you'd like to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Commit your changes.
4. Submit a pull request.


## Acknowledgments

- This project uses the RDKit library for molecular fingerprint generation.
- The Streamlit framework is used for building the user interface.
- The dataset used in this project is derived from public sources.

---
