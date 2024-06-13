import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile, PDBxFile
from pdbfixer.pdbfixer import (
    PDBFixer,
    proteinResidues,
    dnaResidues,
    rnaResidues,
    _guessFileFormat,
)
from flask import (
    Flask,
    request,
    session,
    g,
    render_template,
    make_response,
    send_file,
    redirect,
    url_for,
)
from flask import send_from_directory
import os
from werkzeug.utils import secure_filename
import threading
import time
from pdbfixer.pdbfixer import PDBFixer, _guessFileFormat  # Adjust the import path based on your setup
import sys

if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


from werkzeug.utils import secure_filename
from multiprocessing import Process, Pipe
import datetime
import os
import json
import shutil
import signal
import sys
import tempfile
import threading
import time
import traceback
import webbrowser
import zipfile
import subprocess
import mols2grid
from openff.toolkit import Molecule
from openfe import SmallMoleculeComponent
from openfe.setup import LomapAtomMapper



from flask import send_file
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import io
import base64


# Get the current working directory
current_dir = os.getcwd()


import os
from rdkit.Chem import Draw

uploadedFiles = {}
fixer = None
scriptOutput = None
simulationProcess = None


app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Set a secret key for session management
app.config['UPLOAD_FOLDER'] = 'uploads/'  # Folder to store uploaded files

# Ensure the upload folder exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route("/headerControls")
def headerControls():
    if "startOver" in request.args:
        return showSimulationType()
    if "quit" in request.args:
        func = request.environ.get("werkzeug.server.shutdown")
        if func is None:
            raise RuntimeError("Not running with the Werkzeug Server")
        func()
        return "OpenMM Setup has stopped running.  You can close this window."

@app.route('/')
def index():
    session['fileType'] = 'pdb'  # Set default fileType for demonstration
    return render_template('configureFeFiles.html')

@app.route('/configureFiles', methods=['POST'])
def configureFiles():
    fileType = session.get("fileType", "pdb")
    
    if fileType == "pdb":
        if "file" not in request.files or request.files["file"].filename == "":
            # They didn't select a file. Send them back.
            return render_template('configureFeFiles.html')
        
        # Save the uploaded files
        uploadedFiles = saveUploadedFiles()
        
        configureDefaultOptions()
        
        file, name = uploadedFiles["file"][0]
        file.seek(0, 0)
        session["pdbType"] = _guessFileFormat(file, name)
        
        if session["pdbType"] == "pdb":
            global fixer
            fixer = PDBFixer(pdbfile=file)
        
        preparePDB = request.form.get("preparePDB")
        
        if preparePDB == "no":
            session["sdfFile"] = uploadedFiles["sdfFile"][0][1]  # Save the SDF file path
            return redirect(url_for('display'))
        else:
            return showSelectChains()

    return render_template('configureFeFiles.html')

def saveUploadedFiles():
    uploadedFiles = {}
    for key in request.files:
        file = request.files[key]
        if file.filename != '':
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            if key not in uploadedFiles:
                uploadedFiles[key] = []
            uploadedFiles[key].append((file, filepath))
    return uploadedFiles

def configureDefaultOptions():
    # Placeholder for the actual implementation of configureDefaultOptions
    pass


@app.route('/display')
def display():
    global ligand_mols
    
    # Load ligands using OpenFF toolkit if it's a GET request
    if request.method == 'GET':
        sdf_file_path = session.get("sdfFile")
        
        # Load ligands using OpenFF toolkit
        ligands_sdf = Molecule.from_file(sdf_file_path)

        # Convert these to a list of SmallMoleculeComponent
        ligand_mols = [SmallMoleculeComponent.from_openff(sdf) for sdf in ligands_sdf]
        
        # Generate ligand names for dropdown menu
        ligand_names = [f"{i}_{ligand.name}" for i, ligand in enumerate(ligand_mols)]
        
        im = mols2grid.display(sdf_file_path, size=(200, 200))
        return render_template('display_molecules.html', image=im.data, ligand_names=ligand_names)

    # For POST request, retrieve selected ligands from the form
    elif request.method == 'POST':
        # Retrieve selected ligands from the form
        ligand1_index = int(request.form['ligand1'])
        ligand2_index = int(request.form['ligand2'])

        # Perform further processing with selected ligands

        return redirect(url_for('mapping_page')) 

@app.route('/continue', methods=['POST'])
def continue_simulation():
    return redirect(url_for('mapping_page'))


@app.route('/mapping', methods=['GET', 'POST'])
def mapping_page():
    if request.method == 'POST':
        ligand1_index = int(request.form['ligand1'])
        ligand2_index = int(request.form['ligand2'])

        mapper = LomapAtomMapper()
        lomap_mapping = next(mapper.suggest_mappings(ligand_mols[ligand1_index], ligand_mols[ligand2_index]))

        # Generate and save the figure
        plt.figure(figsize=(10, 5))
        lomap_mapping.display()
        img_io = io.BytesIO()
        plt.savefig(img_io, format='png')
        img_io.seek(0)
        plt.close()

        # Save the image to a file
        img_io.seek(0)
        with open('static/lomap_mapping.png', 'wb') as f:
            f.write(img_io.getbuffer())

        return render_template('mapping.html', ligand_mols=ligand_mols, mapping_image_url='/static/lomap_mapping.png')

    return render_template('mapping.html', ligand_mols=ligand_mols)



@app.route("/selectChains", methods=["POST"])
def selectChains():
    session["heterogens"] = request.form.get("heterogens", "")
    numChains = len(list(fixer.topology.chains()))
    request.form.getlist("include")
    deleteIndices = [
        i for i in range(numChains) if str(i) not in request.form.getlist("include")
    ]
    fixer.removeChains(deleteIndices)
    return showAddResidues()

def showAddResidues():
    spans = []
    chains = list(fixer.topology.chains())
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        return showConvertResidues()
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        chain = chains[key[0]]
        chainResidues = list(chain.residues())
        if key[1] < len(chainResidues):
            offset = int(chainResidues[key[1]].id) - len(residues) - 1
        else:
            offset = int(chainResidues[-1].id)
        spans.append(
            (chain.id, offset + 1, offset + len(residues), ", ".join(residues))
        )
    return render_template("addResidues.html", spans=spans)

@app.route("/addResidues", methods=["POST"])
def addResidues():
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if str(i) not in request.form.getlist("add"):
            del fixer.missingResidues[key]
    return showConvertResidues()

def showConvertResidues():
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        return showAddHeavyAtoms()
    residues = []
    nucleotides = ["DA", "DC", "DG", "DT", "A", "C", "G", "T"]
    for i in range(len(fixer.nonstandardResidues)):
        residue, replaceWith = fixer.nonstandardResidues[i]
        if replaceWith in proteinResidues:
            replacements = proteinResidues
        else:
            replacements = nucleotides
        residues.append(
            (residue.chain.id, residue.name, residue.id, replacements, replaceWith)
        )
    return render_template("convertResidues.html", residues=residues)

@app.route("/convertResidues", methods=["POST"])
def convertResidues():
    for i in range(len(fixer.nonstandardResidues)):
        if str(i) in request.form.getlist("convert"):
            fixer.nonstandardResidues[i] = (
                fixer.nonstandardResidues[i][0],
                request.form["residue" + str(i)],
            )
    fixer.replaceNonstandardResidues()
    return showAddHeavyAtoms()

def showAddHeavyAtoms():
    if session["heterogens"] == "none":
        fixer.removeHeterogens(False)
    elif session["heterogens"] == "water":
        fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    allResidues = list(
        set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys())
    )
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        return addHeavyAtoms()
    residues = []
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        residues.append((residue.chain.id, residue.name, residue.id, ", ".join(atoms)))
    return render_template("addHeavyAtoms.html", residues=residues)

@app.route("/addHeavyAtoms", methods=["POST"])
def addHeavyAtoms():
    fixer.addMissingAtoms()
    return showAddHydrogens()

def showConfigureFiles():
    return render_template('configurePdbFile.html')


@app.route("/getCurrentStructure")
def getCurrentStructure():
    pdb = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb)
    return pdb.getvalue()

def main():

    def open_browser():
        # Give the server a moment to start before opening the browser.
        time.sleep(1)
        url = "http://127.0.0.1:5001"
        webbrowser.open(url)

    threading.Thread(target=open_browser).start()
    app.run(debug=True, port=5001)
    
if __name__ == "__main__":
    main()
