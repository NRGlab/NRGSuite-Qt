# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)
import os.path

#Imports
import argparse
import sys
import pymol
import pandas as pd
from colour import Color

def get_sum_per_residue(surfaces_file):
    residues = []
    values = []
    surf = pd.read_csv(surfaces_file, index_col=0)
    sum_column = surf.sum()
    sum_row = surf.sum(axis='columns')
    for x in range(len(sum_column)):
        values.append(sum_column.iloc[x])
        residues.append(surf.columns[x])
    for y in range(len(sum_row)):
        values.append(sum_row.iloc[y])
        residues.append(surf.index[y])
    return (residues, values)

def get_pairs_contacts(surfaces_file):
    pairs = []
    values = []
    surf = pd.read_csv(surfaces_file, index_col=0)
    for i in range (len(surf.index)):
        for j in range (len(surf.columns)):
            if surf.loc[surf.index[i],surf.columns[j]] != 0:
                pairs.append([surf.index[i],surf.columns[j]])
                values.append(surf.loc[surf.index[i],surf.columns[j]])
    return (pairs, values)

def read_residue(res):
    type_res = res[:3]
    chain_res = res[-1]
    num_res = res[3:-1]
    return (type_res, chain_res, num_res)

def color_residue(res, color,pdb_file):
    type_res, chain_res, num_res = read_residue(res)
    selection_string = os.path.basename(pdb_file[:-4])+'_chain' + chain_res + ' and resi ' + num_res
    pymol.cmd.set_color(res, color)
    pymol.cmd.select('surfaces_sele',selection_string)
    #pymol.cmd.show('spheres', 'surfaces_sele')
    pymol.cmd.set("cartoon_transparency", 0.00, 'surfaces_sele')
    pymol.cmd.color(res, 'surfaces_sele')
    pymol.cmd.delete('surfaces_sele')
    return

def generate_color_scale(values, color_scale_range, color_scale):
    
    if color_scale is None:
        top_color = "red"
        mid_color = "white"
        bottom_color = "blue"
    else:
        color_scale = list(color_scale[1:-1].split(","))
        top_color = color_scale[2]
        mid_color = color_scale[1]
        bottom_color = color_scale[0]
    
    Total_colors = []
    
    for i in range(5):
        c = Color(bottom_color, saturation=1/(i+1))
        Total_colors.append(c.rgb)
    white = Color(mid_color)
    Total_colors.append(white.rgb)
    for i in range(5):
        c = Color(top_color, saturation=1/(5-i))
        Total_colors.append(c.rgb)
    #print (Total_colors)
    
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
        elif n >= len(Total_colors):
            n = len(Total_colors) - 1
        color_codes.append(list(Total_colors[n]))
        
    return (color_codes)

def color_distance(pair, value, color, selected_pairs,pdb_file):
    #create distance object
    distance_string = 'dashed_' + pair[0] + '-' + pair[1]
    type_res1, chain_res1, num_res1 = read_residue(pair[0])
    type_res2, chain_res2, num_res2 = read_residue(pair[1])
    selection_string1 = 'chain ' + chain_res1 + ' and resi ' + num_res1 + ' and n. CA'
    selection_string2 = 'chain ' + chain_res2 + ' and resi ' + num_res2 + ' and n. CA'
    pymol.cmd.set_color(distance_string, color)
    pymol.cmd.distance(distance_string, selection_string1, selection_string2)
    pymol.cmd.group(os.path.basename(pdb_file[:-4])+'_surfaces',distance_string)
    pymol.cmd.color(distance_string, distance_string)
    pymol.cmd.hide('labels', distance_string)
    if pair not in selected_pairs:
        pymol.cmd.disable(distance_string)
    return
    
def label_pairs(pair,selected_pairs,pdb_file):
    #create selection
    pair_string = pair[0] + '-' + pair[1]
    type_res1, chain_res1, num_res1 = read_residue(pair[0])
    type_res2, chain_res2, num_res2 = read_residue(pair[1])
    selection_string1 = os.path.basename(pdb_file[:-4])+'_chain' + chain_res1 + ' & resi ' + num_res1 + ' & n. CA'
    selection_string2 = os.path.basename(pdb_file[:-4])+'_chain' + chain_res2 + ' & resi ' + num_res2 + ' & n. CA'
    #pymol.cmd.select(pair_string, selection_string1 + ' ' + selection_string2)
    #label residues
    pymol.cmd.label(selection_string1,"'%s %s %s' %(resn,resi,chain)")
    pymol.cmd.label(selection_string2,"'%s %s %s' %(resn,resi,chain)")
    selected_residues = pairs_to_residues(selected_pairs)
    if pair[0] not in selected_residues:
        pymol.cmd.hide('labels', selection_string1)
    if pair[1] not in selected_residues:
        pymol.cmd.hide('labels', selection_string2)
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

def split_states(residues, pdb_file):
    chains = []
    for res in residues:
        type_res, chain_res, num_res = read_residue(res)
        if chain_res not in chains:
            chains.append(chain_res)
    for C in chains:
        pymol.cmd.select('surfaces_sele',os.path.basename(pdb_file[:-4])+' and chain ' + C)
        pymol.cmd.extract(os.path.basename(pdb_file[:-4])+'_chain' + C, 'surfaces_sele')
        pymol.cmd.group(os.path.basename(pdb_file[:-4])+'_surfaces',os.path.basename(pdb_file[:-4])+'_chain' + C)
    pymol.cmd.delete(os.path.basename(pdb_file[:-4]))
    return (chains)

def show_separate_surfaces(chains,pdb_file):
    for C in chains:
        pymol.cmd.show('surface', os.path.basename(pdb_file[:-4])+'_chain' + C)
    return
    
def generate_session(pdb_file, surfaces_file, residues_of_interest=None, color_scale=None, color_scale_range=None):
    residues, values = get_sum_per_residue(surfaces_file)
    color_codes = generate_color_scale(values, color_scale_range, color_scale)
    pymol.cmd.load(pdb_file)
    pymol.cmd.color('grey60', os.path.basename(pdb_file[:-4]))
    chains = split_states(residues, pdb_file)
    for C in chains:
        pymol.cmd.set("cartoon_transparency", 0.55, os.path.basename(pdb_file[:-4])+'_chain' + C)
    for i in range(len(residues)):
        if values[i] != 0:
            color_residue(residues[i], color_codes[i],pdb_file)
    pairs, values = get_pairs_contacts(surfaces_file)
    if residues_of_interest is None:
        selected_pairs = get_top_10(pairs, values) # get pairs with largest absolute value of interaction - top 10%
    else:
        residues_of_interest = list(residues_of_interest.split(","))
        selected_pairs = all_pairs_from_interest(pairs, residues_of_interest)
    color_codes = generate_color_scale(values, color_scale_range, color_scale)
    for j in range(len(pairs)):
        color_distance(pairs[j], values[j], color_codes[j], selected_pairs,pdb_file)
    for k in range(len(pairs)):
        label_pairs(pairs[k], selected_pairs,pdb_file)
    show_separate_surfaces(chains,pdb_file)
    pymol.cmd.group('results_surfaces',os.path.basename(pdb_file[:-4])+'_surfaces')
    return


"""
USAGE:
generate_session(pdb_file, input_csv_file, residues_of_interest, color_scale, color_scale_range)
"""
