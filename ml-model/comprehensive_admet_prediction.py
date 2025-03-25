# comprehensive_admet_prediction.py
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
import logging
import datetime
from typing import Dict, Any, List
from fpdf import FPDF
import joblib
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ComprehensiveFeatureExtractor:
    def extract_features(self, smiles_list):
        """Extract comprehensive molecular descriptors from SMILES strings"""
        features = []
        valid_smiles = []
        feature_names = []
        names_captured = False

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                descriptor_dict = {
                    'MolWt': Descriptors.MolWt(mol),
                    'LogP': Descriptors.MolLogP(mol),
                    'TPSA': Descriptors.TPSA(mol),
                    'HBD': Descriptors.NumHDonors(mol),
                    'HBA': Descriptors.NumHAcceptors(mol),
                    'RotBonds': Descriptors.NumRotatableBonds(mol),
                    'AromaticRings': Descriptors.NumAromaticRings(mol),
                    'FractionCSP3': Descriptors.FractionCSP3(mol),
                    'HeavyAtoms': Descriptors.HeavyAtomCount(mol),
                    'Rings': Descriptors.RingCount(mol),
                    'Heteroatoms': Descriptors.NumHeteroatoms(mol),
                    'MolLogP': Descriptors.MolLogP(mol),
                    'MolMR': Descriptors.MolMR(mol),
                    'LabuteASA': Descriptors.LabuteASA(mol),
                    'BalabanJ': Descriptors.BalabanJ(mol),
                    'BertzCT': Descriptors.BertzCT(mol),
                    'Chi0v': Descriptors.Chi0v(mol),
                    'Chi1v': Descriptors.Chi1v(mol),
                    'Chi2v': Descriptors.Chi2v(mol),
                    'Chi3v': Descriptors.Chi3v(mol),
                    'Chi4v': Descriptors.Chi4v(mol),
                    'HallKierAlpha': Descriptors.HallKierAlpha(mol),
                    'Kappa1': Descriptors.Kappa1(mol),
                    'Kappa2': Descriptors.Kappa2(mol),
                    'Kappa3': Descriptors.Kappa3(mol)
                }

                vsa_descriptors = {}
                for i in range(1, 15):
                    vsa_name = f'PEOE_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                for i in range(1, 11):
                    vsa_name = f'SMR_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                for i in range(1, 13):
                    vsa_name = f'SlogP_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                descriptor_dict.update(vsa_descriptors)

                if not names_captured:
                    feature_names = list(descriptor_dict.keys())
                    names_captured = True

                descriptors = [descriptor_dict[name] for name in feature_names]

                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
                fp_bits = [int(b) for b in fp.ToBitString()]

                maccs = MACCSkeys.GenMACCSKeys(mol)
                maccs_bits = [int(b) for b in maccs.ToBitString()]

                features.append(descriptors + fp_bits + maccs_bits)
                valid_smiles.append(smiles)

                if not names_captured:
                    feature_names.extend([f'Morgan_{i}' for i in range(1024)])
                    feature_names.extend([f'MACCS_{i}' for i in range(len(maccs_bits))])
                    names_captured = True

        return np.array(features), valid_smiles, feature_names

class PhysChemAnalyzer:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"Invalid SMILES: {smiles}")

    def calculate_properties(self) -> Dict[str, Any]:
        props = {
            "Molecular Weight": round(Descriptors.MolWt(self.mol), 2),
            "LogP": round(Descriptors.MolLogP(self.mol), 2),
            "H-Bond Donors": Descriptors.NumHDonors(self.mol),
            "H-Bond Acceptors": Descriptors.NumHAcceptors(self.mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(self.mol),
            "TPSA": round(Descriptors.TPSA(self.mol), 2),
            "Aromatic Rings": Descriptors.NumAromaticRings(self.mol),
            "Fraction sp3": round(Descriptors.FractionCSP3(self.mol), 2),
            "QED": round(Descriptors.qed(self.mol), 3)
        }
        props["DrugLikeness"] = "Yes" if (
            props["Molecular Weight"] <= 500 and
            props["LogP"] <= 5 and
            props["H-Bond Donors"] <= 5 and
            props["H-Bond Acceptors"] <= 10
        ) else "No"
        return props

class ComprehensiveADMETPredictor:
    def __init__(self):
        self.models_dir = "admet_models"
        self.feature_extractor = ComprehensiveFeatureExtractor()
        self.models = {}
        self.scalers = {}
        self.model_types = {}
        self._load_models()

    def _load_models(self):
        properties = {
            "caco2": False,
            "pgp_substrate": True,
            "hia": False,
            "bioavailability": False,
            "plasma_binding": False,
            "bbb": True,
            "vdss": False,
            "cyp3a4": True,
            "cyp2d6": True,
            "cyp2c9": True,
            "cyp1a2": True,
            "cyp2c19": True,
            "clearance": False,
            "half_life": False,
            "herg": True,
            "ames": True,
            "hepatotoxicity": True
        }

        for prop, is_classifier in properties.items():
            model_path = os.path.join(self.models_dir, f"{prop}_model.joblib")
            scaler_path = os.path.join(self.models_dir, f"{prop}_scaler.joblib")
            
            if os.path.exists(model_path) and os.path.exists(scaler_path):
                try:
                    self.models[prop] = joblib.load(model_path)
                    self.scalers[prop] = joblib.load(scaler_path)
                    self.model_types[prop] = is_classifier
                except Exception as e:
                    print(f"Error loading {prop} model: {str(e)}")

class ComprehensiveReportGenerator(FPDF):
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=15)
        self.set_margins(10, 10, 10)
        try:
            self.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf', uni=True)
            self.default_font = 'DejaVu'
        except:
            self.default_font = 'Arial'

    def _add_header(self, title):
        self.set_font(self.default_font, 'B', 16)
        self.set_fill_color(20, 40, 80)
        self.set_text_color(255, 255, 255)
        self.cell(0, 10, title, 1, 1, 'C', fill=True)
        self.ln(5)

    def create_report(self, smiles: str, compound_name: str, physchem: Dict, admet: Dict):
        self.add_page()
        self._add_header("Comprehensive Drug Analysis Report")
        
        self.set_font(self.default_font, 'B', 11)
        self.cell(0, 8, f"Compound: {compound_name}", 0, 1, 'C')
        
        self.set_font(self.default_font, '', 9)
        self.cell(0, 6, f"SMILES: {smiles}", 0, 1, 'C')
        self.ln(10)
        
        # Add report content here (similar to original with simplified tables/charts)
        
        report_filename = f"{compound_name.replace(' ', '_')}_report.pdf"
        self.output(report_filename)
        return report_filename

def main():
    print("\n🧪 Comprehensive ADMET Prediction System 🧪\n")
    predictor = ComprehensiveADMETPredictor()

    while True:
        compound_name = input("Enter Compound Name (or 'exit' to quit): ").strip()
        if compound_name.lower() == 'exit':
            break
            
        smiles = input("Enter SMILES string: ").strip()
        if not Chem.MolFromSmiles(smiles):
            print("Invalid SMILES format")
            continue

        try:
            physchem = PhysChemAnalyzer(smiles).calculate_properties()
            features, _, _ = predictor.feature_extractor.extract_features([smiles])
            
            predictions = {}
            for prop, model in predictor.models.items():
                scaled_features = predictor.scalers[prop].transform(features)
                if predictor.model_types[prop]:
                    proba = model.predict_proba(scaled_features)[0, 1]
                    predictions[prop] = "Yes" if proba >= 0.5 else "No"
                else:
                    predictions[prop] = round(model.predict(scaled_features)[0], 2)

            report_gen = ComprehensiveReportGenerator()
            report_file = report_gen.create_report(smiles, compound_name, physchem, predictions)
            print(f"\n✅ Report generated: {report_file}")

        except Exception as e:
            print(f"\n❌ Error: {str(e)}")

if __name__ == "__main__":
    main()