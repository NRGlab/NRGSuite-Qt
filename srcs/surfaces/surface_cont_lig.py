# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

# Imports
import os
import re
import numpy as np

# Useful dicts
aa = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE',
      'N': 'ASN', 'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V': 'VAL', 'E': 'GLU',
      'Y': 'TYR', 'M': 'MET'}


# Input of pdb file and the 2 chains between which we want to evaluate the interactions

def read_residues(pdb_file, chains, ligands):
    list_chains = []
    list_atoms = []
    f = open(pdb_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            if line[21] in chains:
                res_num = re.findall(r'[+-]?\d+', line[22:27])
                res_name = line[17:20]
                string = res_name + str(res_num[0]) + line[21]
                if string not in list_chains and res_name not in ligands:
                    list_chains.append(string)
            if line[17:20] in ligands:
                atom_name = line[12:16].strip()
                atom_number = re.findall(r'[+-]?\d+', line[6:12])
                string = str(atom_number[0]) + atom_name
                if string not in list_atoms:
                    list_atoms.append(string)
    atoms_numbers = []
    for k in range(len(chains)):
        atoms_numbers.append([])
    for line in Lines:
        for i in range(len(chains)):
            if (line[:4] == 'ATOM' or line[:4] == 'HETA') and line[21] == chains[i]:
                atoms_numbers[i].append(int(line[6:12]))
    return (list_chains, list_atoms, atoms_numbers)


# Function to generate the file with the output of perl code vcont.pl

def vcon(pdb_name, custom_vcon_path=None, custom_vcon_out_path=None):
    if custom_vcon_out_path is None:
        vcon_out_path = 'vcon_file.txt'
    else:
        vcon_out_path = custom_vcon_out_path


    if custom_vcon_path is None:
        string = f'{os.path.join(".", "vcon")} {pdb_name} > {custom_vcon_out_path}'
    else:
        string = f'{custom_vcon_path} {pdb_name} > {custom_vcon_out_path}'
    os.system(string)
    return vcon_out_path



# Functions to fix the names of the chains


def get_chain(atom, og_chain, chains, atom_numbers):
    if og_chain == 'L' and atom < 90000:
        for i in range(len(chains)):
            if atom in atom_numbers[i]:
                return (chains[i])
    if og_chain in chains:
        index = chains.index(og_chain)
    else:
        return (og_chain)
    if atom in atom_numbers[index]:
        return (og_chain)
    else:
        for i in range(len(chains)):
            if i != index:
                if atom in atom_numbers[i]:
                    return (chains[i])
    return (0)


# Functions to open vcon results, AMINO.def and MC_st0r5.2_6.dat for atom types interactions

def read_atom(line):
    atnum = int(line[:6])
    attype = line[7:11].strip()
    resnum = int(line[12:17])
    res = line[19:22]
    chain = line[23:24]
    return (atnum, attype, resnum, res, chain)


def read_surface(line):
    surf = (float(line[-12:-6]))
    return (surf)


def read_interactions(file, matrix, chains, ligands, def_file, dat_file, atom_numbers, scale_factor, atoms, res_list):
    f = open(file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if line[:1] != '#' and line != '\n':
            if line[31:34] == 'Sol':
                atnum, attype, resnum, res, chain = read_atom(line)
                fixed_chain = get_chain(atnum, chain, chains, atom_numbers)
                if res in ligands:
                    atom_name = str(atnum) + attype
            else:
                main_atnum, main_attype, main_resnum, main_res, main_chain = read_atom(line[22:])
                fixed_main_chain = get_chain(main_atnum, main_chain, chains, atom_numbers)
                if fixed_main_chain != 0:
                    if main_res not in ligands:
                        main_residue = main_res + str(main_resnum) + fixed_main_chain
                        surf = read_surface(line)
                        if (main_res not in ligands) and (res in ligands):
                            main_residue_index = res_list.index(main_residue)
                            atom_index = atoms.index(atom_name)
                            matrix[main_residue_index, atom_index] += (
                                        surf * score(main_attype, main_res, attype, res, def_file,
                                                     dat_file) * scale_factor)

    return (matrix)


# get the atom type number from def file of choice
def atomtype_num(def_file, res, attyp):
    attyp = attyp.replace(" ", "")
    f = open(def_file, 'r')
    Lines = f.readlines()
    for i in range(len(Lines)):
        if Lines[i][:3] == res:
            ind = Lines[i].index(attyp + ':')
            ind_end = Lines[i][ind:].index(',')
            attype_num = int(Lines[i][ind + len(attyp) + 1:ind + ind_end])
    f.close()
    return (attype_num)


# get the interaction between atom type 1 and atom type 2 from dat file of choice
def interactions(dat_file, at1, at2):
    if len(str(at1)) == 1:
        at1 = ' ' + str(at1)
    if len(str(at2)) == 1:
        at2 = ' ' + str(at2)
    f = open(dat_file, 'r')
    Lines = f.readlines()
    for line in Lines:
        if str(at1) + '-' + str(at2) == line[5:10] or str(at2) + '-' + str(at1) == line[5:10]:
            interact = float(line[13:])
    f.close()
    return (interact)


# get the final score based on atom type num and interactions
def score(attype1, res1, attype2, res2, def_file, dat_file):
    at1 = atomtype_num(def_file, res1, attype1)
    at2 = atomtype_num(def_file, res2, attype2)
    value = interactions(dat_file, at1, at2)
    return (value)


# create file of list of interactions
def list_file(matrix, output_name, atoms, res):
    residues1 = []
    residues2 = []
    values = []
    abs_values = []
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            num = matrix[i][j]
            if num != 0:
                residues1.append(res[i])
                residues2.append(atoms[j])
                values.append(num)
                abs_values.append(abs(num))
    sorted_residues1 = [x for _, x in sorted(zip(abs_values, residues1), reverse=True)]
    sorted_residues2 = [x for _, x in sorted(zip(abs_values, residues2), reverse=True)]
    sorted_values = [x for _, x in sorted(zip(abs_values, values), reverse=True)]
    list_out = os.path.join(os.path.dirname(output_name), 'List_' + os.path.basename(output_name)[:-4] + '.txt')
    f = open(list_out, "w")
    for k in range(len(values)):
        f.write(sorted_residues1[k] + "," + sorted_residues2[k] + "," + str(sorted_values[k]) + "\n")
    f.close()
    return list_out


def read_args():
    import argparse
    parser = argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f", "--pdb_file", action="store")
    parser.add_argument("-c", "--chains", action="store")
    parser.add_argument("-lig", "--ligand", action="store")
    parser.add_argument("-o", "--output_name", action="store")
    parser.add_argument("-def", "--atomtypes_definition", action="store")
    parser.add_argument("-dat", "--atomtypes_interactions", action="store")
    parser.add_argument("-vcon", "--vcon_path", action="store")
    parser.add_argument("-vcon_out", "--vcon_output_path", action="store")
    args = parser.parse_args()
    ligand = args.ligand
    pdb_file = args.pdb_file
    chains = args.chains
    output_name = args.output_name
    vcon_path = args.vcon_path
    vcon_output_path = args.vcon_output_path
    atomtypes_definition = args.atomtypes_definition
    atomtypes_interactions = args.atomtypes_interactions

    main(pdb_file, chains, ligand, output_name, atomtypes_definition, atomtypes_interactions, vcon_path, vcon_output_path)


def fix_csv(output_name, atoms, res):
    with open(output_name) as csv_file:
        lines = csv_file.readlines()
    for line_counter, line in enumerate(lines):
        lines[line_counter] = res[line_counter] + ',' + lines[line_counter]
    atoms_str = ',' + ','.join(atoms) + '\n'
    lines.insert(0, atoms_str)
    return lines


def write_image_file(matrix, output_name, atoms, res):
    out_lines = []
    column = matrix.sum(axis=0)
    row = matrix.sum(axis=1)
    out_lines.append(','.join(atoms) + '\n')
    out_lines.append(','.join(res) + '\n')
    out_lines.append(','.join(np.char.mod('%f', column)) + '\n')
    out_lines.append(','.join(np.char.mod('%f', row)))
    list_out = os.path.join(os.path.dirname(output_name), 'image_' + os.path.basename(output_name)[:-4] + '.txt')
    with open(list_out, "w") as f:
        f.writelines(out_lines)
    return list_out


def main(pdb_file, chains, ligand, output_name, atomtypes_definition, atomtypes_interactions, vcon_path, vcon_output_path):
    list_ligands = ligand.split(",")
    res, atoms, atom_numbers = read_residues(pdb_file, chains, list_ligands)
    # print (res, atoms)
    # print (args.chains, atom_numbers)

    vcon_o_path = vcon(pdb_file, vcon_path, vcon_output_path)

    matrix = np.zeros((len(res), len(atoms)), dtype=np.float64)
    # matrix = np.matrix(temp_matrix)
    # matrix = pd.DataFrame(matrix)
    # matrix.columns = atoms
    # matrix.index = res

    # Determined according to the AB-Bind dataset results
    scale_factor = 0.00024329

    matrix = read_interactions(vcon_o_path, matrix, chains, list_ligands, atomtypes_definition,
                               atomtypes_interactions, atom_numbers, scale_factor, atoms, res)

    # matrix.tofile(output_name, sep=',')
    np.savetxt(output_name, matrix, delimiter=",")
    lines = fix_csv(output_name, atoms, res)
    with open(output_name, 'w') as csv_file:
        csv_file.writelines(lines)
    list_file_path = list_file(matrix, output_name, atoms, res)
    image_file_path = write_image_file(matrix, output_name, atoms, res)
    # remove files
    os.remove(vcon_o_path)

    return list_file_path, image_file_path


if __name__ == '__main__':
    read_args()
