from pymol import cmd
from PyQt5.QtGui import QStandardItemModel, QStandardItem, QColor
import general_functions
import os
from srcs.surfaces import run_Surfaces


def refresh_objects(form):
    model = QStandardItemModel()
    model.setHorizontalHeaderLabels(['Objects'])
    all_objects = cmd.get_object_list('all')
    valid_objects = []
    general_functions.refresh_dropdown_bd_site(form.Settings_wild_type_box,'','',add_none=True, show_all_objects=True)

    # Loop through the lines and add them as items to the model
    for obj in all_objects:
        if cmd.get_type(obj) not in ['group', 'object:measurement', 'selection']:
            valid_objects.append(obj)

    for obj in valid_objects:
        # Create an item with the line text
        item = QStandardItem(obj)
        # Add the item to the model
        model.appendRow(item)
    form.Settings_tableView.setModel(model)
    form.Settings_tableView.setColumnWidth(0, 550)


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
            cmd.group(f'{selected_obj}states',state_obj)


