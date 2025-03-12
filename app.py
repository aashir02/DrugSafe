from flask import Flask, render_template, request, flash
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from sideeffects import SideEffectPredictor
from Ddi_model import predict_side_effects

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Required for flashing messages

# Load SideEffectPredictor model
def load_side_effect_model():
    return SideEffectPredictor.load_model(r"C:\Users\ashir\Music\side_effect_model.pkl")

side_effect_model = load_side_effect_model()

@app.route("/", methods=["GET", "POST"])
def home():
    prediction_results = []  # Default to an empty list
    page = request.form.get("page", "single_drug")  # Default to single drug page

    if request.method == "POST":
        if page == "single_drug":
            smiles_input = request.form.get("smiles_input")
            if smiles_input:
                try:
                    predicted_effects = side_effect_model.predict(smiles_input)
                    if predicted_effects:
                        flash(f"Predicted Side Effects: {predicted_effects}", "success")
                        prediction_results = predicted_effects  # Update prediction_results
                    else:
                        flash("No known side effects found.", "warning")
                except Exception as e:
                    flash(f"Error processing SMILES: {str(e)}", "error")
            else:
                flash("Please enter a valid SMILES structure.", "error")

        elif page == "drug_interaction":
            smiles1 = request.form.get("smiles1")
            smiles2 = request.form.get("smiles2")
            if smiles1 and smiles2:
                try:
                    predicted_interactions = predict_side_effects(smiles1, smiles2)
                    if predicted_interactions:
                        flash(f"Predicted Drug-Drug Interaction Effects: {predicted_interactions}", "success")
                        prediction_results = predicted_interactions  # Update prediction_results
                    else:
                        flash("No known interactions found.", "warning")
                except Exception as e:
                    flash(f"Error processing SMILES structures: {str(e)}", "error")
            else:
                if not smiles1 and not smiles2:
                    flash("Please enter valid SMILES structures for both drugs.", "error")
                elif not smiles1:
                    flash("Please enter a valid SMILES structure for Drug 1.", "error")
                else:
                    flash("Please enter a valid SMILES structure for Drug 2.", "error")

    # Always pass prediction_results to the template
    return render_template("index.html", prediction_results=prediction_results, page=page)

if __name__ == "__main__":
    app.run(debug=True)
