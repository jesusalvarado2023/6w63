import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from stmol import showmol
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
#from padelpy import from_smiles
import numpy as np
#import pickle
import joblib

st.title("Predictor de docking score de ligando-receptor por ML")

compound_smiles=st.text_input('Ingresa tu código SMILES','CC1=CCC(CC1)C(=C)CCC=C(C)C')
mm = Chem.MolFromSmiles(compound_smiles)

Draw.MolToFile(mm,'mol.png')
st.image('mol.png')

#######
RDKit_select_descriptors = joblib.load('./archivos/RDKit_select_descriptors.pickle')
#PaDEL_select_descriptors = joblib.load('./archivos/PaDEL_select_descriptors.pickle')
#robust_scaler = joblib.load('./archivos/robust_scaler.pickle')
#minmax_scaler = joblib.load('./archivos/minmax_scaler.pickle')
#selector_lgbm = joblib.load('./archivos/selector_LGBM.pickle')
#lgbm_model = joblib.load('./archivos/lgbm_best_model.pickle')

# RDKit selected descriptors function
def get_selected_RDKitdescriptors(smile, selected_descriptors, missingVal=None):
    ''' Calculates only the selected descriptors for a molecule '''
    res = {}
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        return {desc: missingVal for desc in selected_descriptors}

    for nm, fn in Descriptors._descList:
        if nm in selected_descriptors:
            try:
                res[nm] = fn(mol)
            except:
                import traceback
                traceback.print_exc()
                res[nm] = missingVal
    return res

df = pd.DataFrame({'smiles': [compound_smiles]})
#st.dataframe(df)

# Calculate selected RDKit descriptors
RDKit_descriptors = [get_selected_RDKitdescriptors(m, RDKit_select_descriptors) for m in df['smiles']]
RDKit_df = pd.DataFrame(RDKit_descriptors)
st.write("Descriptores RDKit")
st.dataframe(RDKit_df)
