from pymol import cmd
import numpy as np
import os
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QPushButton
import shutil
import re

def process_flexaid_result(flexaid_result_file, output):
    with open(flexaid_result_file, 'r') as infile:
        lines = infile.readlines()
    with open(output, 'w') as outfile:
       for line in lines:
            if 'REMARK' not in line:
                if 'LIG  9999' in line:
                    a_name = str(re.sub(r'\d+', '', line[12:17].split()[0])) + str(int(line[9:11])) + ' ' * (
                            5 - len(str(re.sub(r'\d+', '', line[12:17].split()[0])) + str(int(line[9:11]))))
                    new_line = line[:12] + a_name + line[17:21] + 'L' + line[22:]
                    outfile.write(new_line)
                else:
                    outfile.write(line)

def output_message(output_box, text, message_type):
    out_color = None
    red = '<span style="color:red;">{}</span>'
    yellow = '<span style="color:orange;">{}</span>'
    green = '<span style="color:green;">{}</span>'
    if message_type == 'error':
        out_color = red
    elif message_type == 'warning':
        out_color = yellow
    elif message_type == 'valid':
        out_color = green
    output_box.append(out_color.format(text))


def disable_run_mutate_buttons(form, disable=False, enable=False):
    if not disable and not enable:
        print('missing enable or disable args')
        return
    for widget in form.findChildren(QPushButton):
        if "run" in widget.text().lower() or "mutate" in widget.text().lower():
            if disable:
                widget.setDisabled(disable)
            if enable:
                widget.setEnabled(enable)

def show_popup(form, dir_path,temp_path,save_file):

        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        if not save_file:
            msg.setWindowTitle('No sessions found in this project')
            msg.setText('No sessions found in this project,\nDo you want to proceed?')
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        if save_file:
            msg.setWindowTitle('NRGSuite_QT results folder already exists')
            msg.setText('NRGSuite_QT results folder already exists,\nDo you want to proceed?')
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        response = msg.exec_()
        if response == QMessageBox.Ok:
            if not save_file:
                form.temp_line_edit.setText(dir_path)
            if save_file:
                sufix = 2
                base_dir=os.path.join(dir_path, 'NRGSuite_QT_results')
                while os.path.exists(f"{base_dir}_{sufix}"):
                    sufix += 1
                new_dir = f"{base_dir}_{sufix}"
                os.mkdir(new_dir)
                files = os.listdir(temp_path)
                for file_name in files:
                    source_file = os.path.join(temp_path, file_name)
                    target_file = os.path.join(new_dir, file_name)
                    if '.DS_Store' not in file_name:
                        if os.path.isdir(source_file):
                            shutil.copytree(source_file, target_file)
                        else:
                            shutil.copy2(source_file, target_file)
                form.temp_line_edit.setText(os.path.join(new_dir))
                cmd.save(os.path.join(new_dir,'load_project.pse'))
        else:
            show_save_dialog(form, temp_path, save=save_file)

def show_save_dialog(self,temp_path,save=1):

    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    dir_path = QFileDialog.getExistingDirectory(self, "Select Directory", "", options=options)
    if dir_path:
        if save:
            if temp_path!=dir_path:
                files=os.listdir(temp_path)
                if os.path.isdir(os.path.join(dir_path, 'NRGSuite_QT_results')):
                    show_popup(self, dir_path, temp_path, 1)
                else:
                    os.mkdir(os.path.join(dir_path, 'NRGSuite_QT_results'))
                    for file_name in files:
                        source_file = os.path.join(temp_path,file_name)
                        target_file = os.path.join(dir_path,'NRGSuite_QT_results',file_name)
                        if '.DS_Store' not in file_name:
                            if os.path.isdir(source_file):
                                shutil.copytree(source_file, target_file)
                            else:
                                shutil.copy2(source_file, target_file)
                    self.temp_line_edit.setText(os.path.join(dir_path,'NRGSuite_QT_results'))
                    cmd.save(os.path.join(dir_path,'NRGSuite_QT_results','load_project.pse'))
            else:
                cmd.save(os.path.join(dir_path,'load_project.pse'))
        if not save:
            files=os.listdir(dir_path)
            for file in files:
                if file=='load_project.pse':
                    cmd.load(os.path.join(dir_path,'load_project.pse'))
                    self.temp_line_edit.setText(dir_path)
                    break
            else:
                show_popup(self,dir_path,temp_path,0)


def refresh_dropdown_bd_site(dropdown_to_refresh, target, output_box, add_none=False, show_all_objects=False):
    if show_all_objects:
        pymol_objects = sorted(cmd.get_object_list('all'))
        if add_none:
            pymol_objects.insert(0, 'None')
        dropdown_to_refresh.clear()
        dropdown_to_refresh.addItems(pymol_objects)
        dropdown_to_refresh.setCurrentText(pymol_objects[0])
    else:
        if target is None or target == '':
            return
        else:
            loaded_objects = cmd.get_names(type='objects', enabled_only=0)
            if f'gc_{target}' in loaded_objects:
                binding_sites = cmd.get_object_list(f'(gc_{target})')
                if not binding_sites or len(binding_sites) == 0:
                    output_message(output_box, f'No binding sites found for {target}', 'warning')
                else:
                    dropdown_to_refresh.clear()
                    binding_sites = sorted(binding_sites)
                    if add_none:
                        binding_sites.insert(1, 'None')
                    dropdown_to_refresh.addItems(binding_sites)
                    dropdown_to_refresh.setCurrentText(binding_sites[0])
            else:
                output_message(output_box, f'No binding sites found for {target}. You may want to enable "show all objects as binding sites" under settings', 'warning')


def refresh_dropdown_target(dropdown_to_refresh, output_box):
    objects_to_ignore = []
    all_objects = cmd.get_names("all")  # Get all objects and groups in PyMOL
    if len(all_objects) == 0:
        output_message(output_box, f'No object found', 'warning')
    else:
        for pymol_object in all_objects:
            if cmd.get_type(pymol_object) == "object:group":  # Check if it is a group
                objects_to_ignore.append(pymol_object)
                group_members = cmd.get_object_list(f'({pymol_object})')
                for member in group_members:
                    objects_to_ignore.append(member)
        objects_to_display = []
        for pymol_object in all_objects:
            if pymol_object not in objects_to_ignore:
                objects_to_display.append(pymol_object)
        dropdown_to_refresh.clear()
        objects_to_display = sorted(objects_to_display)
        dropdown_to_refresh.addItems(objects_to_display)
        if len(objects_to_display) > 0:
            dropdown_to_refresh.setCurrentText(objects_to_display[0])



def refresh_dropdown(dropdown_to_refresh, output_box, filter_for='', no_warning=False,exclude=None, non_group=1, lig=0, add_none=0):
    list_pymol_objects = cmd.get_names('all')
    if non_group and not lig:
        list_pymol_objects_filtered = cmd.get_object_list('all')
        if 'surfaces_results' in list_pymol_objects:
            list_surfaces=cmd.get_object_list('surfaces_results')
            final_list=[]
            if list_surfaces:
                for obj in list_pymol_objects_filtered:
                    if obj not in list_surfaces:
                        final_list.append(obj)
                list_pymol_objects=final_list
        else:
            list_pymol_objects=list_pymol_objects_filtered
    if lig:
        list_pymol_objects = cmd.get_names('selections')
    if filter_for:
        list_pymol_objects = [x for x in list_pymol_objects if filter_for in x]
    if type(exclude) == list:
        list_pymol_objects = [item for item in list_pymol_objects if all(ex not in item for ex in exclude)]
    if type(exclude) == str:
        list_pymol_objects = [item for item in list_pymol_objects if exclude not in item]
    if len(list_pymol_objects) == 0 and no_warning is False:
        output_message(output_box, 'No objects found', 'warning')
    list_pymol_objects = sorted(list_pymol_objects)
    if add_none:
        list_pymol_objects.insert(1, 'None')
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[0])




def refresh_folder(folder_path, dropdown_to_refresh, ignore_defaults=False):
    # print('===========')
    # print(folder_path)
    folders = next(os.walk(folder_path))[1]
    if ignore_defaults:
        filtered_folders = []
        for folder in folders:
            filenames = next(os.walk(os.path.join(folder_path, folder)), (None, None, []))[2]
            if '.default' not in filenames:
                filtered_folders.append(folder)
        folders = filtered_folders
    folders = sorted([item.replace('_', ' ') for item in folders])
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(folders)


def folder_browser(text_window, ligand_set_path, file_extension):
    smile_file_path = QFileDialog.getOpenFileName(None, 'Select a File', ligand_set_path, file_extension)[0]
    if smile_file_path:
        print(smile_file_path)
        text_window.setText(smile_file_path)


def pymol_hide_structures(form):
    list_pymol_objects = cmd.get_names('all')
    if form.button_hide.isChecked():
        if not list_pymol_objects:
            output_message(form.output_box, 'No clefts to hide', 'warning')
        else:
            form.button_hide.setText('Show')
            cmd.hide('everything', 'bd_site_*')
    else:
        form.button_hide.setText('Hide')
        cmd.show('surface', 'bd_site_*')


def get_mouse_config():
    config_mouse = ''
    try:
        name = cmd.get("button_mode_name")
        if name[0] == '1':
            config_mouse += 'one'
        elif name[0] == '2':
            config_mouse += 'two'
        elif name[0] == '3':
            config_mouse += 'three'
        config_mouse += '_button'
        if name[0] != '1':
            if name[9:] == 'Viewing':
                config_mouse += '_viewing'
            elif name[9:] == 'Editing':
                config_mouse += '_editing'
            elif name[9:] == 'Motions':
                config_mouse += '_motions'
        return config_mouse
    except:
        return 'three_button_viewing'


def read_coords_cleft(cleft_path):
    coords = []
    with open(cleft_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('AT'):
            line = line.split()
            temp_coords = (float(line[6]), float(line[7]), float(line[8]))
            coords.append(temp_coords)
    coords = np.array(coords)
    return lines, coords


def get_residue_info(selection):
    unique_residues = set()
    residue_info = []
    cmd.iterate(selection, 'unique_residues.add((resn, resi, chain))', space={'unique_residues': unique_residues})
    residue_info = [[resname, resn, chain] for resname, resn, chain in unique_residues]
    return residue_info

def surfaces_enable_buttons(form):
    form.flexaid_retrieve_nrgrank_ligands.setEnabled(True)
    form.surfaces_retreive_flexaid_result.setEnabled(True)

def load_color_list(color_list_path):
    with open(color_list_path) as infile:
        lines = infile.readlines()
    list_of_dicts = []
    for line in lines:
        line = line.strip()
        name, rgb = line.split(' ')
        rgb_tuple = tuple(map(int, rgb.strip('()').split(',')))
        list_of_dicts.append({'name': name, 'rgb': rgb_tuple})
    return list_of_dicts


def get_group_of_object(object_name):
    all_objects = cmd.get_names("all")  # Get all objects and groups in PyMOL
    for pymol_object in all_objects:
        if cmd.get_type(pymol_object) == "object:group":  # Check if it is a group
            group_members = cmd.get_object_list(f'({pymol_object})') # Check objects in a group
            if object_name in group_members:
                return pymol_object
    return None