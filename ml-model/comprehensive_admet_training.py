# comprehensive_admet_training.py
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from chembl_webresource_client.new_client import new_client
import joblib
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, mean_absolute_error, roc_auc_score
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

print("\n🧪 Comprehensive ADMET Model Training 🧪\n")

class ComprehensiveFeatureExtractor:
    def extract_features(self, smiles_list):
        """Extract comprehensive molecular descriptors from SMILES strings"""
        features = []
        valid_smiles = []
        feature_names = []

        # Flag to capture feature names once
        names_captured = False

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Calculate comprehensive descriptors
                descriptor_dict = {
                    # Physical properties
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
                    # Lipophilicity and solubility related
                    'MolLogP': Descriptors.MolLogP(mol),
                    'MolMR': Descriptors.MolMR(mol),
                    'LabuteASA': Descriptors.LabuteASA(mol),
                    # Topological
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

                # Add individual VSA descriptors
                vsa_descriptors = {}
                # PEOE_VSA descriptors
                for i in range(1, 15):
                    vsa_name = f'PEOE_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                # SMR_VSA descriptors
                for i in range(1, 11):
                    vsa_name = f'SMR_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                # SlogP_VSA descriptors
                for i in range(1, 13):
                    vsa_name = f'SlogP_VSA{i}'
                    desc_func = getattr(Descriptors, vsa_name, None)
                    if desc_func:
                        vsa_descriptors[vsa_name] = desc_func(mol)

                descriptor_dict.update(vsa_descriptors)

                # Capture feature names once
                if not names_captured:
                    feature_names = list(descriptor_dict.keys())
                    names_captured = True

                # Convert to list in order
                descriptors = [descriptor_dict[name] for name in feature_names]

                # Add Morgan fingerprints
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
                fp_bits = [int(b) for b in fp.ToBitString()]

                # Add MACCS fingerprints
                maccs = MACCSkeys.GenMACCSKeys(mol)
                maccs_bits = [int(b) for b in maccs.ToBitString()]

                # Combine features
                features.append(descriptors + fp_bits + maccs_bits)
                valid_smiles.append(smiles)

                # Extend feature names
                if not names_captured:
                    feature_names.extend([f'Morgan_{i}' for i in range(1024)])
                    feature_names.extend([f'MACCS_{i}' for i in range(len(maccs_bits))])
                    names_captured = True

        return np.array(features), valid_smiles, feature_names

class ComprehensiveDatasetLoader:
    def __init__(self):
        self.data_dir = "admet_datasets"
        os.makedirs(self.data_dir, exist_ok=True)
        self.client = new_client

    def download_chembl_bioactivity(self, target_chembl_id, filename, property_type="IC50", limit=1000):
        """Download bioactivity data from ChEMBL"""
        try:
            activities = self.client.activity.filter(
                target_chembl_id=target_chembl_id,
                standard_type__in=[property_type] if isinstance(property_type, str) else property_type,
                standard_relation__in=["=", "<", ">"],
                standard_value__isnull=False
            ).only(['molecule_chembl_id', 'standard_type', 'standard_value', 'standard_units', 'canonical_smiles'])[:limit]

            df = pd.DataFrame(activities)
            if df.empty:
                return None

            csv_path = os.path.join(self.data_dir, f"{filename}.csv")
            df.to_csv(csv_path, index=False)
            print(f"Downloaded {len(df)} records for {filename}")
            return csv_path
        except Exception as e:
            print(f"Error downloading {filename}: {str(e)}")
            return None

    def download_comprehensive_datasets(self):
        """Download all required datasets"""
        datasets = {}

        # Absorption
        datasets["caco2"] = self.download_chembl_bioactivity("CHEMBL4609", "caco2", "Papp")
        datasets["pgp_substrate"] = self.download_chembl_bioactivity("CHEMBL4302", "pgp_substrate")
        
        # Distribution
        datasets["plasma_binding"] = self.download_chembl_bioactivity("CHEMBL3885882", "plasma_binding", "PPB")
        datasets["bbb"] = self.download_chembl_bioactivity("CHEMBL4523", "bbb", ["logBB", "Kp"])
        
        # Metabolism
        datasets["cyp3a4"] = self.download_chembl_bioactivity("CHEMBL340", "cyp3a4")
        
        # Toxicity
        datasets["herg"] = self.download_chembl_bioactivity("CHEMBL240", "herg")
        
        # Add known drugs
        self._add_known_drugs(datasets)
        return datasets

    def _add_known_drugs(self, datasets):
        """Add known drug data"""
        known_drugs = [
            ["CCOCCN1C(=NC2=CC=CC=C12)C1CCN(CCC2=CC=C(C=C2)C(C)(C)C(O)=O)CC1", "Bilastine", {
                "caco2": 8.5, "pgp_substrate": 1, "plasma_binding": 85, 
                "bbb": 0.01, "cyp3a4": 1, "herg": 0
            }],
            ["CN1C(=NC2=CC=CC=C12)C1CCN(CCC2=CC=C(C=C2)Cl)CC1", "Cetirizine", {
                "caco2": 6.2, "pgp_substrate": 1, "plasma_binding": 93,
                "bbb": 0.06, "cyp3a4": 0, "herg": 0
            }]
        ]

        drugs_dir = os.path.join(self.data_dir, "known_drugs")
        os.makedirs(drugs_dir, exist_ok=True)

        for prop in ["caco2", "pgp_substrate", "plasma_binding", "bbb", "cyp3a4", "herg"]:
            data = []
            for drug in known_drugs:
                if prop in drug[2]:
                    data.append([drug[0], drug[2][prop]])
            
            if data:
                df = pd.DataFrame(data, columns=["smiles", prop])
                path = os.path.join(drugs_dir, f"known_drugs_{prop}.csv")
                df.to_csv(path, index=False)
                datasets[f"known_drugs_{prop}"] = path

class ComprehensiveModelTrainer:
    def __init__(self):
        self.models_dir = "admet_models"
        os.makedirs(self.models_dir, exist_ok=True)
        self.feature_extractor = ComprehensiveFeatureExtractor()

    def train_model(self, datasets, property_name, threshold=None, classification=True):
        """Train individual model"""
        print(f"\nTraining {property_name} model...")
        
        # Prepare data
        X, y, _ = self._prepare_combined_dataset(datasets, property_name, threshold, classification)
        if X is None or len(X) < 10:
            print(f"Insufficient data for {property_name}")
            return None

        # Split and scale
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Model setup
        model = RandomForestClassifier(n_estimators=100, n_jobs=-1) if classification else RandomForestRegressor(n_estimators=100, n_jobs=-1)
        model.fit(X_train, y_train)

        # Evaluate
        if classification:
            proba = model.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, proba)
            print(f"{property_name} AUC: {auc:.3f}")
        else:
            preds = model.predict(X_test)
            mae = mean_absolute_error(y_test, preds)
            print(f"{property_name} MAE: {mae:.3f}")

        # Save
        joblib.dump(model, os.path.join(self.models_dir, f"{property_name}.joblib"))
        joblib.dump(scaler, os.path.join(self.models_dir, f"{property_name}_scaler.joblib"))
        return model

    def _prepare_combined_dataset(self, datasets, property_name, threshold=None, classification=False):
        """Combine data from multiple sources"""
        all_features = []
        all_values = []

        # Main dataset
        if property_name in datasets:
            X, y, _ = self._prepare_dataset(datasets[property_name], 'standard_value', threshold, classification)
            if X is not None:
                all_features.append(X)
                all_values.append(y)

        # Known drugs
        known_key = f"known_drugs_{property_name}"
        if known_key in datasets:
            X, y, _ = self._prepare_dataset(datasets[known_key], property_name, threshold, classification)
            if X is not None:
                all_features.append(X)
                all_values.append(y)

        if not all_features:
            return None, None, None

        return np.vstack(all_features), np.concatenate(all_values), None

    def _prepare_dataset(self, csv_path, value_col, threshold=None, classification=False):
        """Prepare single dataset"""
        try:
            df = pd.read_csv(csv_path)
            smiles = df['canonical_smiles'] if 'canonical_smiles' in df.columns else df['smiles']
            values = df[value_col].values
            
            features, valid_smiles, _ = self.feature_extractor.extract_features(smiles)
            values = values[:len(valid_smiles)]
            
            if classification and threshold is not None:
                values = [1 if float(v) < threshold else 0 for v in values]
            
            return features, np.array(values), valid_smiles
        except Exception as e:
            print(f"Error processing {csv_path}: {str(e)}")
            return None, None, None

if __name__ == "__main__":
    # Initialize components
    loader = ComprehensiveDatasetLoader()
    trainer = ComprehensiveModelTrainer()

    # Download data
    print("\nDownloading datasets...")
    datasets = loader.download_comprehensive_datasets()

    # Train models
    print("\nTraining models...")
    properties = {
        "caco2": False,
        "pgp_substrate": True,
        "plasma_binding": False,
        "bbb": True,
        "cyp3a4": True,
        "herg": True
    }

    for prop, is_classifier in properties.items():
        trainer.train_model(
            datasets=datasets,
            property_name=prop,
            classification=is_classifier,
            threshold=10000 if is_classifier else None
        )

    print("\nTraining complete! Models saved to admet_models/ directory")