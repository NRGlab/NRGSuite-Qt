from pymol import cmd
import numpy as np
from src.getcleft import wizard
from general_functions import read_coords_cleft, output_message
import os


def get_center(cleft_coordinates):
    length = cleft_coordinates.shape[0]
    sum_x = np.sum(cleft_coordinates[:, 0])
    sum_y = np.sum(cleft_coordinates[:, 1])
    sum_z = np.sum(cleft_coordinates[:, 2])
    return [sum_x / length, sum_y / length, sum_z / length]


def get_max_coords(cleft_coordinates, center):
    distances = np.linalg.norm(cleft_coordinates - center, axis=1)
    radius = np.max(distances)
    # max_diameter += 6
    return radius


def display_sphere(cleft_object_name, form, slider, temp_path):
    sphere_name = 'SPHERE'
    if cleft_object_name == '':
        output_message(form.output_box, 'No cleft object selected in step 1', 'warning')
    elif sphere_name in cmd.get_object_list():
        output_message(form.output_box, 'Sphere already exists', 'warning')
    else:
        cmd.hide("everything", cleft_object_name)
        cmd.show("surface", cleft_object_name)
        getcleft_output_path = os.path.join(temp_path, 'GetCleft')
        cleft_save_folder_path = os.path.join(getcleft_output_path, 'Clefts')
        if not os.path.exists(cleft_save_folder_path):
            os.makedirs(cleft_save_folder_path)
        cleft_save_path = os.path.join(cleft_save_folder_path, cleft_object_name + '.pdb')
        cmd.save(cleft_save_path, cleft_object_name)
        cleft_coordinates = np.array(cmd.get_model(cleft_object_name, 1).get_coord_list())
        center_coordinate = get_center(cleft_coordinates)
        max_vdw = get_max_coords(cleft_coordinates, center_coordinate)
        cmd.pseudoatom(sphere_name,
                       pos=center_coordinate,
                       vdw=max_vdw,
                       state=1)
        cmd.color('oxygen', sphere_name)
        cmd.hide('everything', sphere_name)
        cmd.show('spheres', sphere_name)
        cmd.refresh()
        slider_max = int(max_vdw*100)
        slider.setEnabled(True)
        slider.setMaximum(slider_max)
        slider.setValue(slider_max)
        slider.setSingleStep(10)


def resize_sphere(sphere_name, slider_value):
    cmd.alter(sphere_name, 'vdw=' + str(slider_value/100))
    cmd.rebuild(sphere_name)
    cmd.refresh()


def move_sphere():
    cmd.window('hide')
    cmd.window('show')
    wiz = wizard.Sphere()
    cmd.set_wizard(wiz)

def delete_sphere(sphere_object, slider):
    cmd.delete(sphere_object)
    slider.setEnabled(False)

def crop_cleft(sphere_name, sphere_vdw, temp_path, cleft_name, output_box, slider):
    try:
        cmd.count_atoms(sphere_name)
    except Exception:
        output_message(output_box, 'No sphere in PyMOL. Press add in step 2.', 'warning')
        return
    else:
        getcleft_output_path = os.path.join(temp_path, 'GetCleft')
        cleft_save_path = os.path.join(getcleft_output_path, 'Clefts')
        sphere_coords = np.array(cmd.get_model(sphere_name).atom[0].coord)
        min_coords = np.array([sphere_coords-sphere_vdw])
        max_coords = np.array([sphere_coords+sphere_vdw])
        cleft_path = os.path.join(cleft_save_path, cleft_name + '.pdb')
        lines, partition_coords = read_coords_cleft(cleft_path)
        indices = np.where((partition_coords > min_coords).all(axis=1) & (partition_coords < max_coords).all(axis=1))[0] + 1
        lines_to_output = [lines[0]] + [lines[i] for i in indices]
        partition_path = os.path.join(cleft_save_path, cleft_name + '_cropped.pdb')
        with open(partition_path, 'w') as f:
            f.writelines(lines_to_output)
        partition_name = os.path.basename(partition_path).split('.')[0]
        cmd.set("auto_zoom", 0)
        cmd.load(partition_path, format='pdb')
        cmd.hide('everything', partition_name)
        cmd.show('surface', partition_name)
        cmd.color('grey60', partition_name)
        delete_sphere(sphere_name, slider)
        cmd.disable(cleft_name)
        cmd.enable(cleft_name)
        cmd.refresh()

