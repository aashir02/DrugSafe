import streamlit as st
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from sideeffects import SideEffectPredictor
from Ddi_model import predict_side_effects

# Load SideEffectPredictor model
@st.cache_resource
def load_side_effect_model():
    return SideEffectPredictor.load_model(r"C:\Users\ashir\Music\side_effect_model.pkl")

side_effect_model = load_side_effect_model()

# Streamlit App Layout
st.title("💊 Drug Side Effect & Interaction Predictor")
st.sidebar.header("Navigation")
page = st.sidebar.radio("Go to", ["Single Drug Side Effects", "Drug-Drug Interaction"])

if page == "Single Drug Side Effects":
    st.header("Predict Side Effects for a Single Drug")
    smiles_input = st.text_input("Enter a SMILES structure:")
    if st.button("Predict Side Effects"):
        if smiles_input:
            predicted_effects = side_effect_model.predict(smiles_input)
            if predicted_effects:
                st.success("Predicted Side Effects:")
                st.write(predicted_effects)
            else:
                st.warning("No known side effects found.")
        else:
            st.error("Please enter a valid SMILES structure.")

elif page == "Drug-Drug Interaction":
    st.header("Predict Side Effects for a Drug Combination")
    smiles1 = st.text_input("Enter SMILES for Drug 1:")
    smiles2 = st.text_input("Enter SMILES for Drug 2:")
    if st.button("Predict Interactions"):
        if smiles1 and smiles2:
            predicted_interactions = predict_side_effects(smiles1, smiles2)
            if predicted_interactions:
                st.success("Predicted Drug-Drug Interaction Effects:")
                st.write(predicted_interactions)
            else:
                st.warning("No known interactions found.")
        else:
            st.error("Please enter valid SMILES structures for both drugs.")

st.sidebar.markdown("Developed for side effect and drug interaction analysis.")