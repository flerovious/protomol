import os
import pickle
import tempfile

import numpy as np
import streamlit as st
from Bio.PDB import PDBIO, PDBParser

# Load the trained SVM model
model_path = "svm_trained_model.pkl"
with open(model_path, "rb") as file:
    svm_model = pickle.load(file)


def process_pdb(file_path):
    parser = PDBParser()
    structure = parser.get_structure("Protein", file_path)
    modified_structure = structure.copy()
    active_sites = []

    for model in modified_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == "N":
                        coords = np.array([atom.coord])
                        print("Atom coordinates:", coords)
                        prediction = svm_model.predict(coords)
                        print("Prediction:", prediction)
                        if prediction[0] == 1:
                            active_sites.append(residue.id[1])
                            atom.element = "P"  # Modify the element of nitrogen atoms predicted as active
                            print("Active site residue:", residue.id)

    # Save the modified PDB structure
    io = PDBIO()
    io.set_structure(modified_structure)
    temp_dir = tempfile.mkdtemp()
    modified_path = os.path.join(temp_dir, "modified.pdb")
    io.save(modified_path)

    return modified_path, active_sites


st.title("Protein Active Site Identification")

uploaded_file = st.file_uploader("Upload your PDB file", type=["pdb"])
if uploaded_file is not None:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmpfile:
        tmpfile.write(uploaded_file.getvalue())
        modified_path, active_sites = process_pdb(tmpfile.name)

        # Provide download link for the modified PDB file
        with open(modified_path, "rb") as file:
            btn = st.download_button(
                label="Download Modified PDB",
                data=file,
                file_name="modified.pdb",
                mime="application/octet-stream",
            )

        st.write("Active sites residues IDs:", active_sites)
else:
    st.write("Please upload a PDB file to process.")
