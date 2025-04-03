from pymol import cmd
from PyQt5.QtGui import QStandardItemModel, QStandardItem, QColor
import general_functions
import os
from srcs.surfaces import run_Surfaces


def refresh_hetatm(form, silence_warnings=False):
    selection = form.Settings_tableView.selectionModel().selectedRows()
    if len(selection) != 1:
        general_functions.output_message(form.output_box, 'Only one object may be selected at a time to use this function', 'warning')
    else:
        selected_obj = form.Settings_tableView.model().data(selection[0])
        dropdown_to_refresh = form.Settings_list_hetatm
        hetatms = sorted(list(set(atom.resn for atom in cmd.get_model(f"{selected_obj} and hetatm").atom)))
        if len(hetatms) > 0:
            hetatms.append('All')
            dropdown_to_refresh.clear()
            dropdown_to_refresh.addItems(hetatms)
            dropdown_to_refresh.setCurrentIndex(0)
        else:
            if not silence_warnings:
                general_functions.output_message(form.output_box, 'No objects found', 'warning')


def remove_het_atoms(form):
    selection = form.Settings_tableView.selectionModel().selectedRows()
    if len(selection) != 1:
        general_functions.output_message(form.output_box,
                                         'Only one object may be selected at a time to use this function', 'warning')
    else:
        selected_obj = form.Settings_tableView.model().data(selection[0])
        molecule_to_remove = form.Settings_list_hetatm.currentText( )
        if molecule_to_remove == 'All':
            cmd.remove(f"hetatm and {selected_obj}")
        else:
            cmd.remove(f'resn {molecule_to_remove} and {selected_obj}')
        refresh_hetatm(form, silence_warnings=True)


############

def refresh_objects(form):
    all_objects = cmd.get_object_list('all')
    if len(all_objects) == 0:
        general_functions.output_message(form.output_box, 'No objects found', 'warning')
    else:
        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(['Objects'])
        valid_objects = []
        general_functions.refresh_dropdown_bd_site(form.Settings_wild_type_box,'','',add_none=True, show_all_objects=True)

        for obj in all_objects:
            if cmd.get_type(obj) not in ['group', 'object:measurement', 'selection']:
                valid_objects.append(obj)

        for obj in valid_objects:
            item = QStandardItem(obj)
            model.appendRow(item)
        form.Settings_tableView.setModel(model)
        form.Settings_tableView.setColumnWidth(0, 550)
        form.tool_group_1.setEnabled(True)
        form.tool_group_2.setEnabled(True)
        form.tool_group_3.setEnabled(True)


def combine_objects(form):
    selection = form.Settings_tableView.selectionModel().selectedRows()
    selected_obj=[]
    for index in selection:
        selected_obj.append(form.Settings_tableView.model().data(index))
    if form.Settings_combine_obj_name.text() == '':
        general_functions.output_message(form.output_box, 'Please define a name for the multi_state obj', 'error')
    else:
        for obj in selected_obj:
                cmd.join_states(form.Settings_combine_obj_name.text(), obj)

def devide_states(form):
    selection = form.Settings_tableView.selectionModel().selectedRows()
    selected_obj=[]
    if form.Settings_wild_type_box.currentText() != 'None':
        target = form.Settings_wild_type_box.currentText()
        target_file = os.path.join(form.obj_manager_save_path, target + '.pdb')
        cmd.save(target_file, target)
    for index in selection:
        selected_obj.append(form.Settings_tableView.model().data(index))
    for obj in selected_obj:
        cmd.split_states(obj)
        for state in range(cmd.count_states(obj)):
            state_obj = obj + "_{:04}".format(state+1)
            if form.Settings_wild_type_box.currentText() != 'None':
                    state_obj=obj + "_{:04}".format(state+1)
                    state_file=os.path.join(form.obj_manager_save_path, obj + "_{:04}.pdb".format(state+1))
                    cmd.save(state_file,state_obj)
                    diff=run_Surfaces.compare_residues(target_file,state_file)
                    cmd.set_name(obj+"_{:04}".format(state+1),obj+f'_{diff}')
                    state_obj=obj+f'_{diff}'
                    os.remove(state_file)
            cmd.group(f'{obj}states',state_obj)
