import concurrent.futures
import os
import sys
install_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(install_dir)
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from src.nrgrank.nrgrank_general_functions import get_params_dict
from itertools import repeat
from datetime import datetime
from rdkit.Chem import rdDistGeom
from src.nrgrank.process_ligands import preprocess_ligands_one_target
import subprocess
from pathlib import Path
import json


def load_rad_dict(filepath):
    with open(filepath, 'r') as file:
        loaded_atom_data = json.load(file)
    return loaded_atom_data


def generate_conformers(smiles_line, no_conformers, name_position):
    etkdg = rdDistGeom.ETKDGv3()
    # === optional settings ===
    # etkdg.maxAttempts = 10
    # etkdg.pruneRmsThresh = 0.5
    # etkdg.numThreads = 10
    # https://greglandrum.github.io/rdkit-blog/posts/2023-03-02-clustering-conformers.html
    etkdg.randomSeed = 0xa700f
    etkdg.verbose = False
    etkdg.useRandomCoords = True  # Start with random coordinates
    split_smiles_line = smiles_line.split()
    smiles = split_smiles_line[0]
    name = split_smiles_line[name_position]
    molecule = Chem.MolFromSmiles(smiles)
    try:
        frags = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=False)
    except:
        print('Error getting fragment for: ', smiles_line)
        frags = molecule
        if frags is None:
            return None
    molecule = max(frags, key=lambda frag: frag.GetNumAtoms())
    # print(smiles_line)
    mol_weight = rdMolDescriptors.CalcExactMolWt(molecule)
    num_heavy_atoms = molecule.GetNumHeavyAtoms()
    # print(f"Molecular weight of the molecule in Daltons: {mol_weight:.2f} Da")
    if mol_weight > 1800 or num_heavy_atoms <= 3:
        print("Molecular weight is over 1800 Da or the molecule has fewer than 3 heavy atoms: Ignoring")
        return None
    else:
        mol = Chem.AddHs(molecule, addCoords=True)
        if no_conformers == 1:
            try:
                AllChem.EmbedMolecule(mol, params=etkdg)
            except Exception as e:
                print('=====================================')
                print(f'Error: {e}\n Molecule: {smiles_line}\n')
                print('=====================================')
                return None
        else:
            AllChem.EmbedMultipleConfs(mol, no_conformers, params=etkdg)
        mol.SetProp("_Name", name)
    return mol


def read_params():
    smiles_path = sys.argv[1]
    optimize = sys.argv[2]
    custom_output_file_path = sys.argv[3]
    try:
        custom_deps_path = sys.argv[4]
    except IndexError:
        custom_deps_path = None
    convert = True
    preprocess = True
    main(smiles_path, optimize, custom_output_file_path, preprocess, convert, custom_deps_path)


def main(smiles_path, optimize, custom_output_file_path, preprocess, convert, custom_deps_path):
    print("Started generating conformers @ ", datetime.now())
    root_software_path = Path(__file__).resolve().parents[1]
    os.chdir(root_software_path)
    name_position = -1
    print("Dont forget to specify in what column the name is for the smiles file")
    print("Path to smiles: ", smiles_path)
    config_path = os.path.join(root_software_path, "deps/config.txt")
    if custom_deps_path:
        config_path = os.path.join(custom_deps_path, "config.txt")
    params_dict = get_params_dict(config_path)
    conf_num = params_dict["CONFORMER_NUMBER"]
    if custom_output_file_path != "False":
        sdf_output_file = custom_output_file_path
        output_folder = os.path.dirname(sdf_output_file)
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
    else:
        output_folder = os.path.join(os.path.dirname(smiles_path), f"{os.path.basename(smiles_path).split('.')[0]}_conformers")
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        end = "_not_opt"
        if optimize == "yes":
            end = "_opt"
        sdf_output_file = os.path.join(output_folder, f"{os.path.splitext(os.path.basename(smiles_path))[0]}_{conf_num}_conf{end}.sdf")
    mol2_output_file = os.path.splitext(sdf_output_file)[0] + '.mol2'

    writer = AllChem.SDWriter(sdf_output_file)
    if conf_num == 0:
        exit("number of conformers is 0")
    with open(smiles_path) as f:
        lines = f.readlines()
        if lines[0].startswith('smiles'):
            lines = lines[1:]

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for mol in executor.map(generate_conformers, lines, repeat(conf_num), repeat(name_position)):
            if mol is not None:
                for cid in range(mol.GetNumConformers()):
                    if optimize == "yes":
                        Chem.rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol, cid)
                    mol = Chem.RemoveHs(mol)
                    writer.write(mol, cid)

    print("done with writter")
    AllChem.SDWriter.close(writer)
    print("Finished generating conformers @ ", datetime.now())
    if convert:
        print("converting to mol2")
        open_babel_command = f"obabel \"{sdf_output_file}\" -O \"{mol2_output_file}\" ---errorlevel 1"
        print(f'obabel command: {open_babel_command}')
        subprocess.run(open_babel_command, shell=True, check=True)
        print("removing")
        os.remove(sdf_output_file)
    if preprocess:
        rad_dict_path = os.path.join(root_software_path, "deps", "atom_type_radius.json")
        if custom_deps_path:
            rad_dict_path = os.path.join(custom_deps_path, "atom_type_radius.json")
        rad_dict = load_rad_dict(rad_dict_path)
        preprocess_ligands_one_target(rad_dict, conf_num, output_folder,
                                      'single_file', mol2_output_file,
                                      None)
        # process_ligands.main("single_file", None, output_folder, mol2_output_file, None)


if __name__ == '__main__':
    read_params()
