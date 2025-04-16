import os
from nrgten.encom import ENCoM
from nrgten.atom import Atom
import argparse
import json
import numpy as np
import sys
import pickle
import time
install_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(install_dir)
from srcs.surfaces.ligand_atomtypes import add_pdb
from srcs.surfaces.run_Surfaces import create_ligand_file, flex_res
from general_functions import process_flexaid_result
# vcon compiled on Windows with the developper command prompt breaks vcon


def generate_massfile(pdb_filename, mass_filename):
        atoms = []
        with open(pdb_filename) as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM') or line.startswith("HETATM"):
                atoms.append(Atom(line))
        xyz_data = np.zeros((len(atoms), 3))
        for i, atom in enumerate(atoms):
            xyz_data[i] = atom.xyz
        centroid = np.array([np.mean(xyz_data[:, 0]), np.mean(xyz_data[:, 1]), np.mean(xyz_data[:, 2])])
        medoid = None
        mindist = float('Inf')
        other_atom = None
        other_flag = True
        for atom in atoms:
            dist = np.linalg.norm(atom.xyz - centroid)
            if dist < mindist:
                mindist = dist
                medoid = atom
            elif other_flag:
                other_atom = atom
                other_flag = False
        if other_flag:
            other_atom = atoms[0]
        with open(mass_filename, "w") as f:
            f.write("CONNECT: {} -> {}\n".format(medoid.name, other_atom.name))
            f.write("N_MASSES: 1\n")
            f.write("CENTER: {}\n".format(medoid.name))
            f.write("NAME: {}\n".format(medoid.name))


def find_het( target_file, temp_path, main_folder_path):
    het_dic = {}
    with open(target_file, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            if 'HETATM' in line:
                het_dic[line[17:20]] = '1'
    for lig in list(het_dic):
        def_file = os.path.join(main_folder_path, "deps", "surfaces", 'AMINO_FlexAID.def')
        open_def_file = open(def_file, "r")
        ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
        create_ligand_file(target_file, ligand_file_name)
        custom_def_path = os.path.join(temp_path, 'NRGTEN', f'custom_{os.path.basename(def_file)}')
        custom_def_file = open(custom_def_path, 'w')
        add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
        open_def_file.close()
        custom_def_file.close()
        generate_massfile(ligand_file_name + '.pdb', os.path.join(temp_path, 'NRGTEN', f'{lig}.masses'))
        with open(custom_def_path, "r") as t2:
            texto = t2.readlines()
            definition = texto[-1][:-2] + '\n'
        with open(os.path.join(temp_path, 'NRGTEN', f'{lig}.atomtypes'), "w") as t3:
            t3.write(definition)
    return list(het_dic)


def flex_aid_matrix(main_folder_path):
    matrix = np.array([np.zeros(41)] * 41)
    flexaid_dat_path = os.path.join(main_folder_path, "deps", "surfaces", "FlexAID.dat")
    with open(flexaid_dat_path, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            matrix[int(line.split("=")[0].split('-')[0])][int(line.split("=")[0].split('-')[1])] = -float(
                line.split("=")[1])
    return matrix


def encom_model(target_file, main_folder_path, temp_path, list_het=None):
    if flex_res(target_file):
        process_flexaid_result(target_file, target_file)
    if not list_het:
        list_het = find_het(target_file, temp_path, main_folder_path)
    atype_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.atomtypes')]
    mass_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.masses')]
    for lig in list_het:
        atype_lists.append(os.path.join(temp_path, 'NRGTEN', f'{lig}.atomtypes'))
        mass_lists.append(os.path.join(temp_path, 'NRGTEN', f'{lig}.masses'))

    matrix = flex_aid_matrix(main_folder_path)
    model = ENCoM(target_file, interact_mat=matrix, atypes_list=atype_lists, massdef_list=mass_lists)
    return model


def write_b_factor(target, dyna_sig, temp_path, labels):
    target_file = os.path.join(temp_path, 'NRGTEN', f'{target}.pdb')
    b_factor_dict = {}
    with open(os.path.join(temp_path, 'NRGTEN', target + '_dynasig.txt'), 'w') as t1:
        for res in range(len(labels)):
            t1.write('{} {}\n'.format(labels[res], dyna_sig[res]))
    for i in range(len(dyna_sig)):
        key = '{}_{}_{}'.format(labels[i].split('|')[0][:3], labels[i].split('|')[2], labels[i].split('|')[1])
        b_factor_dict[key] = dyna_sig[i]
    with open(target_file, 'r') as t1:
        texto = t1.readlines()
        texto_list = []
        for line in texto:
            if 'HETATM' in line or 'ATOM' in line:
                key = '{}_{}_{}'.format(line[17:20], line[21], int(line[22:26]))
                number = b_factor_dict[key]
                lined = f"{number:.2f}"
                lined_abs = f"{abs(number):.2f}"
                lined = lined.rjust(6)
                lined_abs = lined_abs.rjust(6)
                texto_list.append(line[:54] + lined + lined_abs + line[66:])
            else:
                texto_list.append(line)
    output_path = os.path.join(temp_path, 'NRGTEN', f'{target}_dynasig.pdb')
    with open(output_path, 'w') as t2:
        for line in texto_list:
            t2.write(line)
    return b_factor_dict


def run_dynamical_signature(target_file, beta, main_folder_path, temp_path, list_het=None):
    start_time = time.time()
    target = os.path.splitext(os.path.basename(target_file))[0]
    model = encom_model(target_file, main_folder_path, temp_path, list_het)
    dyna_sig = model.compute_bfactors_boltzmann(beta=float(beta))
    svib = model.compute_vib_entropy(beta=float(beta))
    model_mass_label = model.get_mass_labels()
    b_fact_dictionary = write_b_factor(target, dyna_sig, temp_path, model_mass_label)
    pickle_file_path = os.path.splitext(target_file)[0] + '.pkl'
    with open(pickle_file_path, "wb") as f:
        pickle.dump((b_fact_dictionary, dyna_sig, model_mass_label, svib), f)
    end_time = time.time()
    execution_time = end_time - start_time
    with open(os.path.splitext(pickle_file_path)[0] + '.txt' , "w") as f:
        f.write(f"Execution time: {execution_time:.4f} seconds")
    print(f"Execution time: {execution_time:.4f} seconds")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--target_file', required=True, type=str, help='Path to target_file')
    parser.add_argument('-b', '--beta', type=str, help='beta')
    parser.add_argument('-m', '--main_folder_path', type=str, help='main_folder_path')
    parser.add_argument("-te", '--temp_path', type=str, help='Path to temp')
    parser.add_argument("-l", '--list_het', type=str, help='list_het')

    args = parser.parse_args()

    target_file = args.target_file
    beta = args.beta
    main_folder_path = args.main_folder_path
    temp_path = args.temp_path
    list_het = args.list_het
    list_het= json.loads(list_het)
    run_dynamical_signature(target_file, beta, main_folder_path, temp_path, list_het=list_het)


if __name__ == '__main__':
    get_args()