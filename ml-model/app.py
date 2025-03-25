from flask import Flask, request, jsonify
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import os  # <-- Add this line!

app = Flask(__name__)

# Example ML prediction function (replace with your actual model)
def predict_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "num_atoms": mol.GetNumAtoms()
    }

@app.route('/predict', methods=['POST'])
def predict():
    data = request.json
    smiles = data.get('smiles', '')
    results = predict_properties(smiles)
    if not results:
        return jsonify({"error": "Invalid SMILES string"}), 400
    return jsonify(results)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 8080))  # Now works with `os` imported
    app.run(host='0.0.0.0', port=port)