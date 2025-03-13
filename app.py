from flask import Flask, render_template, request, flash, redirect, url_for, session
import pandas as pd
import pickle
import json
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from sideeffects import SideEffectPredictor
from Ddi_model import predict_side_effects
from functools import wraps

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Required for flashing messages and sessions

# User database file
USERS_FILE = 'users.json'

# Load SideEffectPredictor model
def load_side_effect_model():
    return SideEffectPredictor.load_model(r"C:\Users\ashir\Music\side_effect_model.pkl")

side_effect_model = load_side_effect_model()

# Create users.json if it doesn't exist
def init_users_file():
    if not os.path.exists(USERS_FILE):
        with open(USERS_FILE, 'w') as f:
            json.dump({}, f)

# Load users from JSON file
def load_users():
    init_users_file()
    try:
        with open(USERS_FILE, 'r') as f:
            return json.load(f)
    except:
        return {}

# Save users to JSON file
def save_users(users):
    with open(USERS_FILE, 'w') as f:
        json.dump(users, f)

# Login required decorator
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'username' not in session:
            flash('Please log in to access this page', 'error')
            return redirect(url_for('login'))
        return f(*args, **kwargs)
    return decorated_function

@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        username = request.form.get("username")
        password = request.form.get("password")
        
        users = load_users()
        
        if username in users and users[username] == password:
            session['username'] = username
            flash(f"Welcome back, {username}!", "success")
            return redirect(url_for('home'))
        else:
            flash("Invalid username or password", "error")
    
    return render_template("login.html")

@app.route("/register", methods=["GET", "POST"])
def register():
    if request.method == "POST":
        username = request.form.get("username")
        password = request.form.get("password")
        confirm_password = request.form.get("confirm_password")
        
        if not username or not password:
            flash("Username and password are required", "error")
        elif password != confirm_password:
            flash("Passwords don't match", "error")
        else:
            users = load_users()
            if username in users:
                flash("Username already exists", "error")
            else:
                users[username] = password
                save_users(users)
                flash("Registration successful! Please log in.", "success")
                return redirect(url_for('login'))
    
    return render_template("register.html")

@app.route("/logout")
def logout():
    session.pop('username', None)
    flash("You have been logged out", "success")
    return redirect(url_for('login'))
# Add these imports to the top of your app.py
from flask import send_file
import csv
import io
from datetime import datetime

# Add this new route to your app.py
@app.route("/download_report", methods=["POST"])
@login_required
def download_report():
    # Get the prediction results from the form
    prediction_results = request.form.getlist('prediction_results[]')
    drug_type = request.form.get('drug_type', 'drug')
    smiles_input = request.form.get('smiles_input', '')
    smiles1 = request.form.get('smiles1', '')
    smiles2 = request.form.get('smiles2', '')
    
    # Create an in-memory file for the CSV
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Write the header
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    username = session.get('username', 'Unknown User')
    
    if drug_type == 'single_drug':
        writer.writerow(['Drug Side Effect Report'])
        writer.writerow(['Generated on:', timestamp])
        writer.writerow(['User:', username])
        writer.writerow(['SMILES Structure:', smiles_input])
        writer.writerow([])
        writer.writerow(['Predicted Side Effects:'])
    else:  # drug_interaction
        writer.writerow(['Drug Interaction Report'])
        writer.writerow(['Generated on:', timestamp])
        writer.writerow(['User:', username])
        writer.writerow(['Drug 1 SMILES:', smiles1])
        writer.writerow(['Drug 2 SMILES:', smiles2])
        writer.writerow([])
        writer.writerow(['Predicted Interaction Effects:'])
    
    # Write the results
    for result in prediction_results:
        writer.writerow([result])
    
    # Reset the pointer to the beginning of the file
    output.seek(0)
    
    # Generate a filename
    if drug_type == 'single_drug':
        filename = f"drug_side_effects_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    else:
        filename = f"drug_interaction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    
    # Return the file as an attachment
    return send_file(
        io.BytesIO(output.getvalue().encode('utf-8')),
        mimetype='text/csv',
        as_attachment=True,
        download_name=filename
    )

@app.route("/", methods=["GET", "POST"])
@login_required
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

    # Get the username from the session
    username = session.get('username')
    
    # Always pass prediction_results and username to the template
    return render_template("index.html", prediction_results=prediction_results, page=page, username=username)

if __name__ == "__main__":
    app.run(debug=True)
