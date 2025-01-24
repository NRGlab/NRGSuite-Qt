# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)


import argparse
import os.path

# ANY ATOM NOT IN THIS LIST WILL BE CONSIDERED DUMMY:39
dict_types={'C.1':'1', 'C.2':'2', 'C.3':'3', 'C.ar':'4', 'C.cat':'5', 'N.1':'6',
 'N.2':'7', 'N.3':'8', 'N.4':'9', 'N.ar':'10', 'N.am':'11', 'N.pl3':'12',
 'O.2':'13', 'O.3':'14', 'O.co2':'15', 'O.ar':'16', 'S.2':'17', 'S.3':'18',
 'S.o':'19', 'S.O2':'20', 'S.ar':'21', 'P.3':'22', 'F':'23', 'Cl':'24', 'Br':'25', 'I':'26', 'Se':'27', 'Mg':'28', 'Sr':'29',
 'Cu':'30', 'Mn':'31', 'Hg':'32', 'Cd':'33', 'Ni':'34', 'Zn':'35', 'Ca':'36', 'Fe':'37', 'Cooh':'38', 'Du':'39'}

import os
import pandas as pd
import argparse
import pymol

def convertto_mol2(pdb_file):
    pymol.cmd.load(pdb_file)
    # Updated to use os.path.splitext() and os.path.join()
    mol2_file = os.path.splitext(pdb_file)[0] + '.mol2'
    pymol.cmd.save(mol2_file, os.path.basename(os.path.splitext(pdb_file)[0]))
    pymol.cmd.delete(os.path.basename(os.path.splitext(pdb_file)[0]))
    return

def convertto_pdb(mol_file):
    pymol.cmd.load(mol_file)
    # Updated to use os.path.splitext()
    pdb_file = os.path.splitext(mol_file)[0] + '.pdb'
    pymol.cmd.save(pdb_file)
    return

def read_mol2(mol_file):
    with open(mol_file, 'r') as f:
        Lines = f.readlines()
    selected_Lines = []
    list_atomnames = []
    list_atomtypes = []
    list_res = []
    select = False
    for i in range(len(Lines)):
        if Lines[i] == '@<TRIPOS>BOND\n':
            select = False
        if select:
            selected_Lines.append(Lines[i])
        if Lines[i] == '@<TRIPOS>ATOM\n':
            select = True
    for line in selected_Lines:
        add = line[:-2].split("\t")
        print(add)
        list_atomnames.append(add[1])
        list_atomtypes.append(add[-4])
        list_res.append(add[-2])
    for k in range(len(list_atomnames)):
        list_atomnames[k] = list_atomnames[k].strip()
    res = list_res[0]
    print('res:', res)
    return (list_atomnames, list_atomtypes, res)

def atomtypes_to_numbers(list_atomtypes):
    list_atomnumbers = []
    for atomtype in list_atomtypes:
        if atomtype in dict_types:
            list_atomnumbers.append(dict_types[atomtype])
        else:
            list_atomnumbers.append(dict_types['Du'])
    return list_atomnumbers

def check_atoms(list_atomnames, list_atomtypes):
    for i in range(len(list_atomnames)):
        for j in range(len(list_atomnames)):
            if list_atomnames[i] == list_atomnames[j] and list_atomtypes[i] != list_atomtypes[j]:
                return False
    return True

def check_atom_names_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        Lines = f.readlines()
    for line in Lines:
        if (line[:4] == 'ATOM' or line[:4] == 'HETA'):
            atom_name = line[12:17].strip()
            if len(atom_name) > 3:
                return False
    return True

def check_atom_names_mol2(mol2_file):
    with open(mol2_file, 'r') as f:
        Lines = f.readlines()
    selected_Lines = []
    select = False
    for i in range(len(Lines)):
        if Lines[i] == '@<TRIPOS>BOND\n':
            select = False
        if select:
            selected_Lines.append(Lines[i])
        if Lines[i] == '@<TRIPOS>ATOM\n':
            select = True
    df = pd.DataFrame(columns=['number','name','coord1','coord2','coord3','type','n','res','k'])
    for line in selected_Lines:
        add = line[:-2].split("\t")
        df.loc[len(df)] = add
    list_atomnames = df['name'].tolist()
    for atom_name in list_atomnames:
        atom_name = atom_name.strip()
        if len(atom_name) > 3:
            return False
    return True

def custom_def_file(f, g, list_atomnames, list_atomnumbers, res):
    Lines = f.readlines()
    new_line = '\n' + res[:3] + ' | '
    for i in range(len(list_atomnames)):
        new_line += list_atomnames[i] + ':' + list_atomnumbers[i] + ', '
    for line in Lines:
        g.write(line)
    g.write(new_line)
    return

def remove_duplicates(list_atomnames, list_atomnumbers):
    new_list_atomnames = []
    new_list_atomnumbers = []
    for k in range(len(list_atomnames)):
        if list_atomnames[k] not in new_list_atomnames:
            new_list_atomnames.append(list_atomnames[k])
            new_list_atomnumbers.append(list_atomnumbers[k])
    return new_list_atomnames, new_list_atomnumbers

def add_mol2(ligand_file, f, g):
    list_atomnames, list_atomtypes, res = read_mol2(ligand_file)
    if not check_atom_names_mol2(ligand_file):
        print("WARNING: ATOM NAMES LARGER THAN 3 CHARACTERS - POSSIBLE PROBLEM WITH ATOM TYPE READING")
    if not check_atoms(list_atomnames, list_atomtypes):
        print("WARNING: ATOMS WITH DIFFERENT ATOM TYPES AND SAME ATOM NAME")
    list_atomnumbers = atomtypes_to_numbers(list_atomtypes)
    list_atomnames, list_atomnumbers = remove_duplicates(list_atomnames, list_atomnumbers)
    custom_def_file(f, g, list_atomnames, list_atomnumbers, res)
    return

def add_pdb(ligand_file, f, g):
    convertto_mol2(ligand_file)
    mol2_file = os.path.splitext(ligand_file)[0] + '.mol2'
    list_atomnames, list_atomtypes, res = read_mol2(mol2_file)
    if not check_atom_names_pdb(ligand_file):
        print("WARNING: ATOM NAMES LARGER THAN 3 CHARACTERS - POSSIBLE PROBLEM WITH ATOM TYPE READING")
    if not check_atoms(list_atomnames, list_atomtypes):
        print("WARNING: ATOMS WITH DIFFERENT ATOM TYPES AND SAME ATOM NAME")
    list_atomnumbers = atomtypes_to_numbers(list_atomtypes)
    list_atomnames, list_atomnumbers = remove_duplicates(list_atomnames, list_atomnumbers)
    custom_def_file(f, g, list_atomnames, list_atomnumbers, res)
    return

def main():
    parser = argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-pdb", "--ligand_pdb_file", action="store")
    parser.add_argument("-mol2", "--ligand_mol2_file", action="store")
    parser.add_argument("-def", "--atomtypes_definition", action="store")
    parser.add_argument("-prefix", "--prefix_new_atomtypes_definition", action="store")
    args = parser.parse_args()

    f = open(args.atomtypes_definition, 'r')
    if args.prefix_new_atomtypes_definition is None:
        g = open('custom_' + args.atomtypes_definition, 'w')
    else:
        g = open(args.prefix_new_atomtypes_definition + '_' + args.atomtypes_definition, 'w')

    if args.ligand_mol2_file is None:
        list_ligand_files = args.ligand_pdb_file.split(",")
        for ligand_file in list_ligand_files:
            add_pdb(ligand_file, f, g)
    else:
        list_ligand_files = args.ligand_mol2_file.split(",")
        for ligand_file in list_ligand_files:
            add_mol2(ligand_file, f, g)

    g.close()
    return

if __name__ == '__main__':
    main()
