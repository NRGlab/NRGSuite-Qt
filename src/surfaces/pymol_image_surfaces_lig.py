# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

#Imports
import argparse
import sys
import pymol
import os


def read_image_data(output_name):
    with open(output_name, "r") as f:
        lines = f.readlines()
    atoms = lines[0].strip().split(',')
    residues = lines[1].strip().split(',')
    values_atoms = list(map(float, lines[2].strip().split(',')))
    value_residues = list(map(float, lines[3].strip().split(',')))
    return residues, atoms, value_residues, values_atoms


def get_pairs_contacts(surfaces_file):
    pairs = []
    values = []
    with open(surfaces_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip().split(',')
        pairs.append([line[0], line[1]])
        values.append(float(line[2]))
    return pairs, values

def read_residue(res):
    type_res = res[:3]
    chain_res = res[-1]
    num_res = res[3:-1]
    return (type_res, chain_res, num_res)

def read_atom(atom):
    atom_num = ''
    atom_name = ''
    num = True
    for i in range(len(atom)):
        if atom[i].isnumeric() and num:
            atom_num = atom_num + atom[i]
        else:
            num = False
            atom_name = atom_name + atom[i]
    return (atom_name, atom_num)

def color_residue(res, color, pdb_file):
    type_res, chain_res, num_res = read_residue(res)
    selection_string = os.path.basename(pdb_file[:-4])+'_chain' + chain_res + ' and resi ' + num_res
    pymol.cmd.set_color(res, color)
    pymol.cmd.select('sele_surf',selection_string)
    #pymol.cmd.show('spheres', 'sele')
    pymol.cmd.set("cartoon_transparency", 0.00, 'sele_surf')
    pymol.cmd.color(res, 'sele_surf')
    pymol.cmd.delete('sele_surf')
    return


def load_total_colors(file_path):
    tc = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip().split(',')
        tc.append((float(line[0]), float(line[1]), float(line[2])))
    return tc


def generate_color_scale(values, color_scale_range, color_scale, color_rgb_file_path):
    
    if color_scale is None:
        top_color = "red"
        mid_color = "white"
        bottom_color = "blue"
    else:
        color_scale = list(color_scale[1:-1].split(","))
        top_color = color_scale[2]
        mid_color = color_scale[1]
        bottom_color = color_scale[0]
    # TODO: add ability to generate custum color file
    Total_colors = load_total_colors(color_rgb_file_path)

    
    if color_scale_range is None:
        max_value = max(values)
        min_value = min(values)
        if abs(min_value) > abs(max_value):
            range_value = 2 * abs(min_value)
        else:
            range_value = 2 * abs(max_value)
        step_value = range_value/10
    else:
        min_value = float(color_scale_range[0])
        max_value = float(color_scale_range[1])
        range_value = max_value - min_value
        step_value = range_value/10
    
    color_codes = []
    for value in values:
        s = range_value/2 - (-1*value)
        n = int(s // step_value)
        if n < 0:
            n = 0
        elif n > len(Total_colors):
            n = len(Total_colors)
        color_codes.append(list(Total_colors[n]))
        
    return (color_codes)

def color_distance(pair, value, color, selected_pairs,pdb_file):
    #create distance object
    distance_string = os.path.basename(pdb_file[:-4])+'_dashed_' + pair[0] + '-' + pair[1]
    distance_string = distance_string.replace("'", "")
    type_res, chain_res, num_res = read_residue(pair[0])
    atom_name, atom_num = read_atom(pair[1])
    selection_string1 = os.path.basename(pdb_file[:-4])+'_chain' + chain_res + ' and resi ' + num_res + ' and n. CA'
    selection_string2 = 'id ' + atom_num +' and '+os.path.basename(pdb_file[:-4])+'_ligand'
    pymol.cmd.set_color(distance_string, color)
    pymol.cmd.distance(distance_string, selection_string1, selection_string2)
    pymol.cmd.color(distance_string, distance_string)
    pymol.cmd.group(os.path.basename(pdb_file[:-4]) + '_surfaces',distance_string)
    pymol.cmd.hide('labels', distance_string)
    if pair not in selected_pairs:
        pymol.cmd.disable(distance_string)
    return
    
def label_pairs(pair,selected_pairs,pdb_file):
    #create selection
    pair_string = pair[0] + '-' + pair[1]
    type_res, chain_res, num_res = read_residue(pair[0])
    atom_name, atom_num = read_atom(pair[1])
    selection_string1 = os.path.basename(pdb_file[:-4])+'_chain' + chain_res + ' and resi ' + num_res + ' and n. CA'
    selection_string2 = 'id ' + atom_num +' and '+os.path.basename(pdb_file[:-4])+'_ligand'
    pymol.cmd.select(pair_string, selection_string1 + ' ' + selection_string2)
    #label residues
    pymol.cmd.label(selection_string1,"'%s %s %s' %(resn,resi,chain)")
    selected_residues = pairs_to_residues(selected_pairs)
    if pair[0] not in selected_residues:
        pymol.cmd.hide('labels', selection_string1)
    pymol.cmd.disable(pair_string)
    pymol.cmd.set("label_position", "[0,0,6]")
    pymol.cmd.delete(pair_string)
    return
    
def pairs_to_residues(pairs):
    residues = []
    for i in range(len(pairs)):
        for j in range(len(pairs[0])):
            if pairs[i][j] not in residues:
                residues.append(pairs[i][j])
    return (residues)

def get_top_10(pairs, values):
    top_pairs = []
    absolute_values = []
    size_10_percent = len(values)//10
    for value in values:
        absolute_values.append(abs(value))
    absolute_values.sort(reverse=True)
    top_values = absolute_values[:size_10_percent]
    for f in range(len(pairs)):
        if len(top_pairs) <= len(top_values):
            if (values[f] in top_values) or (-1*values[f] in top_values):
                top_pairs.append(pairs[f])
    return (top_pairs)

def all_pairs_from_interest(pairs, residues_of_interest):
    selected_pairs = []
    for pair in pairs:
        if pair[0] in residues_of_interest or pair[1] in residues_of_interest:
            selected_pairs.append(pair)
    return (selected_pairs)

def split_states(residues, atoms, pdb_file):
    chains = []
    atom_nums = []
    for res in residues:
        type_res, chain_res, num_res = read_residue(res)
        if chain_res not in chains:
            chains.append(chain_res)
    for atom in atoms:
        atom_name, atom_num = read_atom(atom)
        atom_nums.append(atom_num)
    selection_string = ''
    for num in atom_nums:
        selection_string = selection_string + num + '+'
    pymol.cmd.select('sele_surfaces','id ' + selection_string[:-1] + ' and '+os.path.basename(pdb_file)[:-4])
    pymol.cmd.extract(os.path.basename(pdb_file[:-4])+'_ligand', 'sele_surfaces')
    pymol.cmd.group(os.path.basename(pdb_file[:-4])+'_surfaces',os.path.basename(pdb_file[:-4])+'_ligand')
    for C in chains:
        pymol.cmd.select('sele_surfaces','chain ' + C +' and '+os.path.basename(pdb_file[:-4]))
        pymol.cmd.extract(os.path.basename(pdb_file[:-4])+'_chain' + C, 'sele_surfaces')
        pymol.cmd.group(os.path.basename(pdb_file[:-4]) + '_surfaces',os.path.basename(pdb_file[:-4])+'_chain' + C)
    pymol.cmd.delete(os.path.basename(pdb_file[:-4]))
    pymol.cmd.delete('sele_surfaces')
    return (chains)

def show_separate_surfaces(chains,pdb_file):
    for C in chains:
        pymol.cmd.show('surface', os.path.basename(pdb_file[:-4])+'_chain' + C)
        #pymol.cmd.set('transparency', 0.7, 'chain' + C)
    return

def color_ligands(pdb_file):
    pymol.cmd.color("cyan",os.path.basename(pdb_file[:-4])+'_ligand')
    pymol.util.cnc(os.path.basename(pdb_file[:-4])+'_ligand')
    return


def generate_session(pdb_file, image_file, list_file, color_rgb_file_path, residues_of_interest=None, color_scale=None, color_scale_range=None):
    residues, atoms, values_residues, values_atoms = read_image_data(image_file)
    color_codes = generate_color_scale(values_residues, color_scale_range, color_scale, color_rgb_file_path)
    pymol.cmd.load(pdb_file)
    pymol.cmd.color('grey60', os.path.basename(pdb_file)[:-4])
    chains = split_states(residues, atoms, pdb_file)
    for C in chains:
         pymol.cmd.set("cartoon_transparency", 0.55, os.path.basename(pdb_file[:-4])+'_chain' + C)
    for i in range(len(residues)):
         if values_residues[i] != 0:
             color_residue(residues[i], color_codes[i],pdb_file)
    pairs, values = get_pairs_contacts(list_file)
    if residues_of_interest is None:
        selected_pairs = pairs
    else:
        residues_of_interest = list(residues_of_interest.split(","))
        selected_pairs = all_pairs_from_interest(pairs, residues_of_interest)
    color_codes = generate_color_scale(values, color_scale_range, color_scale, color_rgb_file_path)
    for j in range(len(pairs)):
        color_distance(pairs[j], values[j], color_codes[j], selected_pairs, pdb_file)
    for k in range(len(pairs)):
        label_pairs(pairs[k], selected_pairs,pdb_file)
    pymol.cmd.group('surfaces_results',os.path.basename(pdb_file[:-4]) + '_surfaces')
    #show_separate_surfaces(chains,pdb_file)
    color_ligands(pdb_file)
    return

"""
USAGE:
generate_session(pdb_file, input_csv_file, residues_of_interest, color_scale, color_scale_range)
"""
