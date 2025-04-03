from general_functions import get_residue_info
import os
import matplotlib.pyplot as plt
from srcs.surfaces.ligand_atomtypes import add_pdb
from general_functions import process_flexaid_result
from srcs.surfaces.clean_structure import main as clean_structure
from srcs.surfaces.surface_cont_lig import main as surface_cont_lig
from srcs.surfaces.surface_cont import main as surface_cont
from srcs.surfaces.pymol_image_surfaces_lig import generate_session
from srcs.surfaces.pymol_image_surfaces import generate_session as generate_session_ppi
from pymol import cmd
import pandas as pd
from PyQt5.QtGui import QStandardItemModel, QStandardItem, QColor
import csv
from Bio import PDB
from matplotlib.colors import LinearSegmentedColormap
import time
import general_functions

def refresh_res(form,out_path):
    l_dir=os.listdir(out_path)
    ind_res = [file[:-4] for file in l_dir if 'List_' in file]
    cf_comp= [file[:-4] for file in l_dir if '_diff.csv' in file]
    form.surface_select_individual_result.clear()
    form.surface_select_individual_result.addItems(ind_res)
    load_csv_data(form,os.path.join(out_path,form.surface_select_individual_result.currentText()+'.txt'))
    form.surface_select_cf_comparison.clear()
    form.surface_select_cf_comparison.addItems(cf_comp)


def load_csv_data(form, csv_file):

    model = QStandardItemModel()

    with open(csv_file, newline='', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        data = list(csv_reader)
    if len(data[0]) == 2:
        headers = ['MODEL', 'CF']
        last_column_index = 1
    else:
        headers = ['RESIDUE', 'RESIDUE', 'CF']
        last_column_index = 2

    model.setHorizontalHeaderLabels(headers)

    if len(data) > 0:
        last_column_values = [float(row[last_column_index]) for row in data]
        min_value = min(last_column_values)
        max_value = max(last_column_values)

        normalize = lambda x: 0.5 + (x / (2 * max(abs(min_value), abs(max_value)))) if max_value > min_value else 0.5

        colors = [(0, 0, 1), (0, 0, 0), (1, 0, 0)]  # Blue -> Black -> Red
        cmap = LinearSegmentedColormap.from_list("custom_bkr", colors)
        for row_idx, row in enumerate(data):
            items = []
            row[-1] = str(round(float(row[-1]), 4))
            for col_idx, field in enumerate(row):
                item = QStandardItem(field)
                if col_idx == last_column_index:
                    normalized_value = normalize(float(field))
                    color = QColor(*[int(c * 255) for c in cmap(normalized_value)[:3]])
                    item.setBackground(color)
                items.append(item)
            model.appendRow(items)

    form.surfaces_tableView.setModel(model)

    if len(headers)>0:
        form.surfaces_tableView.setColumnWidth(last_column_index, 175)
    form.surfaces_tableView.clicked.connect(lambda index: open_res(form.surfaces_tableView, index,len(headers)))


def get_residues_from_pdb(file_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    residues = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip HETATMs or water molecules
                if PDB.is_aa(residue, standard=True):
                    res_id = residue.get_id()[1]
                    chain_id = chain.get_id()
                    res_name = residue.get_resname()
                    residues[(chain_id, res_id)] = res_name
    return residues


def compare_residues(pdb1, pdb2):
    residues_pdb1 = get_residues_from_pdb(pdb1)
    residues_pdb2 = get_residues_from_pdb(pdb2)

    differences = []

    for res_id, res_name_pdb1 in residues_pdb1.items():
        res_name_pdb2 = residues_pdb2.get(res_id)
        if res_name_pdb2 and res_name_pdb1 != res_name_pdb2:
            chain, res_num = res_id
            differences.append(f"{res_name_pdb1}{res_num}{res_name_pdb2}_{chain}")

    return '_'.join(differences)


def open_res(tableView, index, lheaders):
    if lheaders==3:
        if index.column()==0 or index.column()==1:
            cell_text = tableView.model().data(index)
            cmd.select('sele_surfaces','resname {} and resi {} and chain {}'.format(cell_text[:3],cell_text[3:-1],cell_text[-1]))
            cmd.zoom('sele_surfaces', buffer=10.0)
            cmd.show('lines', 'sele_surfaces')


def create_ligand_file(pdb_file_name, lig_path):
    lig_name = os.path.basename(lig_path).strip()
    with open(pdb_file_name, "r") as f:
        lines = f.readlines()
    lig_pdb_file = open(lig_path + ".pdb", "w")
    for line in lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            res = line[17:20].strip()
            if res == lig_name:
                lig_pdb_file.write(line)
        if 'CONECT' in line:
            lig_pdb_file.write(line)
    lig_pdb_file.close()
    return


def flex_res(target_file):
    with open(target_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if 'LIG  9999' in line:
                return 1
    return 0



def get_chains_from_object(object_name):
    chains = cmd.get_chains(object_name)
    str_chain = ''
    for chain in chains:
        str_chain = str_chain + str(chain)
    return str_chain


def load_surfaces(form, temp_path, main_folder_path, binary_folder_path, binary_suffix):
    general_functions.output_message(form.output_box, '=========== Surfaces ===========', 'valid')
    start_time = time.time()
    vcon_binary_path = os.path.join(binary_folder_path, f'vcon{binary_suffix}')
    target = form.surface_select_object_1.currentText()
    if target == '':
        general_functions.output_message(form.output_box, 'No object selected', 'warning')
        return
    lig = form.surface_select_ligand_object_1.currentText()
    target_2 = form.surface_select_object_2.currentText()
    lig_2 = form.surface_select_ligand_object_2.currentText()
    chain_1=form.surfaces_chain.text()
    chain_2=form.surfaces_chain_2.text()

    if lig != 'None':
        if target_2 == 'None':
            target_file =  os.path.join(temp_path,'Surfaces', f'{target}_{lig}.pdb')
            cmd.save(target_file, target)
            ligands = get_residue_info(lig)
            for ligand in ligands:
                target_chain = get_chains_from_object(target)
                interac_dic = run_surfaces_lig(target_file, target_chain, ligand[0], temp_path, main_folder_path, vcon_binary_path,form)
        else:
            cf_dic = {}
            target_file =  os.path.join(temp_path,'Surfaces', f'{target}_{lig}.pdb')
            cmd.save(target_file, target)
            ligands = get_residue_info(lig)
            for ligand in ligands:
                target_chain = get_chains_from_object(target)
                interac_dic = run_surfaces_lig(target_file, target_chain, ligand[0], temp_path, main_folder_path, vcon_binary_path,form)
                interac_dic = dict(sorted(interac_dic.items(), key=lambda item: abs(item[1]), reverse=True))
                cf_ref = sum(interac_dic.values())
            for state in range(cmd.count_states(target_2)):
                output_file = os.path.join(temp_path, f'{target_2}_{state + 1}.pdb')
                cmd.save(output_file, target_2, state=state + 1)
                diff = compare_residues(target_file, output_file) + f'_{chain_1}_{chain_2}'
                os.rename(output_file, os.path.join(temp_path, 'Surfaces', f'{target_2}_{diff}.pdb'))
                output_file = os.path.join(temp_path, 'Surfaces', f'{target_2}_{diff}.pdb')
                target_chain = get_chains_from_object(target_2)
                ligands = get_residue_info(lig_2)
                for ligand in ligands:
                    target_chain = get_chains_from_object(target)
                    interac_dic = run_surfaces_lig(output_file, target_chain, ligand[0], temp_path, main_folder_path, vcon_binary_path,form)
                    cf_dic[os.path.basename(output_file) + '_' + ligand[0]] = sum(interac_dic.values()) - cf_ref
            csv_file=os.path.join(temp_path,'Surfaces','{}_{}'.format(target,target_2)+'_diff.csv')
            with open(csv_file,'w') as t1:
                for item in list(cf_dic):
                    t1.write('{},{}\n'.format(item,cf_dic[item]))
            load_csv_data(form,csv_file)

    else:
        if chain_1!='None' and chain_2!='None':
            if target_2 == 'None':
                target_file =  os.path.join(temp_path,'Surfaces', f'{target}_{chain_1}_{chain_2}.pdb')
                cmd.save(target_file, target)
                csv_file_surf=run_surfaces_ppi(target_file, chain_1,chain_2,temp_path, main_folder_path, vcon_binary_path,form)
                read_and_select_residues(csv_file_surf,target,num_rows=10)
            else:
                target_file =  os.path.join(temp_path,'Surfaces', f'{target}_{chain_1}_{chain_2}.pdb')
                cmd.save(target_file, target)
                csv_file_surf = run_surfaces_ppi(target_file, chain_1, chain_2, temp_path, main_folder_path,
                                                 vcon_binary_path, form)
                #read_and_select_residues(csv_file_surf, target, num_rows=10)
                cf_ref=cf_calculatior(csv_file_surf)
                cf_dic={}
                for state in range(cmd.count_states(target_2)):
                    output_file = os.path.join(temp_path,'Surfaces', f'{target_2}_{chain_1}_{chain_2}_{state + 1}.pdb')
                    cmd.save(output_file, target_2, state=state + 1)
                    diff = compare_residues(target_file, output_file) + f'_{chain_1}_{chain_2}'
                    os.rename(output_file, os.path.join(temp_path, 'Surfaces', f'{target_2}_{diff}.pdb'))
                    output_file = os.path.join(temp_path, 'Surfaces', f'{target_2}_{diff}.pdb')
                    csv_file_surf = run_surfaces_ppi(output_file, chain_1, chain_2, temp_path, main_folder_path,
                                                 vcon_binary_path, form)
                    cf_dic[os.path.basename(output_file)]=cf_calculatior(csv_file_surf)-cf_ref
                    os.rename(csv_file_surf, os.path.join(temp_path,'Surfaces', f'List_{target}_{diff}.txt'))
                csv_file = os.path.join(temp_path, 'Surfaces', '{}_{}'.format(target, target_2) + '_diff.csv')
                with open(csv_file, 'w') as t1:
                    for item in list(cf_dic):
                        t1.write('{},{}\n'.format(item, cf_dic[item]))
                read_and_select_residues(csv_file, target)
    cmd.disable(f'GetCleft,NRGRank,FlexAID')
    end_time = time.time()
    execution_time = end_time - start_time
    general_functions.output_message(form.output_box, f"Execution time: {execution_time:.4f} seconds", 'valid')
    general_functions.output_message(form.output_box, '=========== END Surfaces ===========', 'valid')
    form.Surfaces_tabs.setCurrentIndex(1)


def cf_calculatior(res_file):
    with open(res_file, 'r') as t1:
        texto=t1.readlines()
        interaction=0
        for line in texto:
            interaction+=float(line.split(',')[-1])
    return interaction


def read_and_select_residues(file_path, object_name, select_first_column_only=True, num_rows=None):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    selection_list = []

    if num_rows is not None:
        if num_rows != 'ALL':
            num_rows=int(num_rows)
            lines = lines[:num_rows]
            print(f'TOP {num_rows} selected in all_residues')

    for line in lines:

        parts = line.strip().split(',')
        if len(parts) >= 2:
            residue1 = parts[0]
            residue2 = parts[1]

            res1_name, res1_id, res1_chain = residue1[:3], residue1[3:-1], residue1[-1]

            selection_criteria1 = f"({object_name} and resn {res1_name} and resi {res1_id} and chain {res1_chain})"

            selection_list.append(selection_criteria1)

            if select_first_column_only:
                res2_name, res2_id, res2_chain = residue2[:3], residue2[3:-1], residue2[-1]

                selection_criteria2 = f"({object_name} and resn {res2_name} and resi {res2_id} and chain {res2_chain})"

                selection_list.append(selection_criteria2)

    combined_selection = " or ".join(selection_list)
    cmd.select('all_residues', combined_selection)


def run_surfaces_ppi(target_file,chain_1, chain_2, temp_path, main_folder_path, vcon_binary_path,form):
    surfaces_output_path = os.path.join(temp_path, 'Surfaces')
    def_file = os.path.join(main_folder_path, "deps", "surfaces", 'AMINO_FlexAID.def')
    flexaid_dat_path = os.path.join(main_folder_path, "deps", "surfaces", 'FlexAID.dat')
    if flex_res(target_file):
        process_flexaid_result(target_file, target_file)
    cleaned_file_path = os.path.join(os.path.dirname(target_file), f"cleaned_{os.path.basename(target_file)}")
    clean_pdb_file = open(cleaned_file_path, "w")
    clean_structure(target_file, def_file, clean_pdb_file)
    clean_pdb_file.close()
    os.rename(cleaned_file_path,
              os.path.join(os.path.dirname(target_file), f'{os.path.splitext(os.path.basename(target_file))[0]}_.pdb'))
    cleaned_file_path=os.path.join(os.path.dirname(target_file), f'{os.path.splitext(os.path.basename(target_file))[0]}_.pdb')
    output_name=os.path.join(surfaces_output_path,os.path.basename(target_file)[:-4]+'.csv')
    surface_cont(cleaned_file_path,chain_1,chain_2,output_name,def_file,flexaid_dat_path,vcon_binary_path)
    load_csv_data(form,os.path.join(os.path.dirname(output_name),"List_" + os.path.basename(output_name[:-4]) + ".txt"))
    generate_session_ppi(cleaned_file_path,output_name)
    return os.path.join(os.path.dirname(output_name),"List_" + os.path.basename(output_name[:-4]) + ".txt")


def run_surfaces_lig(target_file, target_chain, lig, temp_path, main_folder_path, vcon_path,form):
    surfaces_output_path = os.path.join(temp_path, 'Surfaces')
    def_file = os.path.join(main_folder_path, "deps", "surfaces", 'AMINO_FlexAID.def')
    flexaid_dat_path = os.path.join(main_folder_path, "deps", "surfaces", 'FlexAID.dat')
    color_rgb_path = os.path.join(main_folder_path, "deps", "surfaces", 'color_rgb.txt')
    open_def_file = open(def_file, "r")
    if flex_res(target_file):
        process_flexaid_result(target_file, target_file)
    ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
    create_ligand_file(target_file, ligand_file_name)
    custom_def_path = os.path.join(surfaces_output_path, f'custom_{os.path.basename(def_file)}')
    custom_def_file = open(custom_def_path, 'w')
    add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
    open_def_file.close()
    custom_def_file.close()
    cleaned_file_path = os.path.join(os.path.dirname(target_file), f"cleaned_{os.path.basename(target_file)}")
    clean_pdb_file = open(cleaned_file_path, "w")
    clean_structure(target_file, custom_def_path, clean_pdb_file)
    clean_pdb_file.close()
    os.rename(cleaned_file_path, os.path.join(os.path.dirname(target_file), f'{os.path.splitext(os.path.basename(target_file))[0]}_.pdb'))
    cleaned_file_path = os.path.join(os.path.dirname(target_file), f'{os.path.splitext(os.path.basename(target_file))[0]}_.pdb')
    vcon_out_file = os.path.join(os.path.dirname(surfaces_output_path), 'vcon_file.txt')
    csv_path = os.path.join(surfaces_output_path, f'{os.path.basename(target_file)[:-4]}_csv_output.csv')
    list_file_path, image_file_path = surface_cont_lig(cleaned_file_path, target_chain, lig, csv_path, custom_def_path, flexaid_dat_path, vcon_path, vcon_out_file)
    load_csv_data(form, list_file_path)
    generate_session(cleaned_file_path, image_file_path, list_file_path, color_rgb_path)
    interact_dic = {}
    with open(list_file_path, 'r') as f:
        texto = f.readlines()
        for line in texto:
            if line.split(',')[0][3:-1] in interact_dic.keys():
                interact_dic[line.split(',')[0][3:-1]] += float(line.split(',')[-1][:-1])
            else:
                interact_dic[line.split(',')[0]] = float(line.split(',')[-1][:-1])
    return interact_dic


def plot_interactive_table(labels, values):
    df = pd.DataFrame({'Res': labels, 'Value': values})
    fig, ax = plt.subplots()
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.auto_set_column_width(col=[0, 1])
    plt.show()




