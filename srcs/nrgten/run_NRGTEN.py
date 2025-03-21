from nrgten.encom import ENCoM
from nrgten.atom import Atom
from pymol import cmd
from srcs.surfaces.ligand_atomtypes import add_pdb
from srcs.surfaces.run_Surfaces import create_ligand_file, flex_res, compare_residues
from general_functions import process_flexaid_result
import os
import numpy as np
import plotly.graph_objects as go
import time
from general_functions import output_message


def remove_selection_and_save(object_name, selection, output_file):
    cmd.create("object_without_selection", f"{object_name} and not ({selection})")
    cmd.save(output_file, "object_without_selection")
    cmd.delete("object_without_selection")


def flex_aid_matrix(main_folder_path):
    matrix = np.array([np.zeros(41)] * 41)
    flexaid_dat_path =  os.path.join(main_folder_path, "deps","surfaces", "FlexAID.dat")
    with open(flexaid_dat_path, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            matrix[int(line.split("=")[0].split('-')[0])][int(line.split("=")[0].split('-')[1])] = -float(
                line.split("=")[1])
    return matrix


def encom_model(target_file, main_folder_path, temp_path):
    if flex_res(target_file):
        process_flexaid_result(target_file, target_file)
    list_het = find_het(target_file, temp_path, main_folder_path)
    atype_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.atomtypes')]
    mass_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.masses')]
    for lig in list_het:
        atype_lists.append(os.path.join(temp_path, 'NRGTEN', f'{lig}.atomtypes'))
        mass_lists.append(os.path.join(temp_path,'NRGTEN', f'{lig}.masses'))

    matrix = flex_aid_matrix(main_folder_path)
    model = ENCoM(target_file, interact_mat=matrix, atypes_list=atype_lists, massdef_list=mass_lists)
    return model


def standardize_to_minus1_plus1(data):
    max_abs_value = max(abs(x) for x in data)
    standardized_data = [x / max_abs_value for x in data]
    return standardized_data


def find_het(target_file, temp_path, main_folder_path):
    het_dic = {}
    with open(target_file, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            if 'HETATM' in line:
                het_dic[line[17:20]] = '1'
    for lig in list(het_dic):
        def_file = os.path.join(main_folder_path, "deps","surfaces", 'AMINO_FlexAID.def')
        flexaid_dat_path = os.path.join(main_folder_path, "deps","surfaces", 'FlexAID.dat')
        open_def_file = open(def_file, "r")
        ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
        create_ligand_file(target_file, ligand_file_name)
        custom_def_path = os.path.join(temp_path,'NRGTEN', f'custom_{os.path.basename(def_file)}')
        custom_def_file = open(custom_def_path, 'w')
        add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
        open_def_file.close()
        custom_def_file.close()
        generate_massfile(ligand_file_name + '.pdb', os.path.join(temp_path, 'NRGTEN', f'{lig}.masses'))
        with open(custom_def_path, "r") as t2:
            texto = t2.readlines()
            definition = texto[-1][:-2] + '\n'
        with open(os.path.join(temp_path,'NRGTEN', f'{lig}.atomtypes'), "w") as t3:
            t3.write(definition)
    return list(het_dic)


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


def write_b_factor(target, dyna_sig, temp_path, labels):
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    b_factor_dict = {}
    with open(os.path.join(temp_path, 'NRGTEN',target+'_dynasig.txt'),'w') as t1:
        for res in range(len(labels)):
            t1.write('{} {}\n'.format(labels[res],dyna_sig[res]))
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
                lined_abs=f"{abs(number):.2f}"
                lined = lined.rjust(6)
                lined_abs= lined_abs.rjust(6)
                texto_list.append(line[:54] + lined + lined_abs+ line[66:])
            else:
                texto_list.append(line)
    output_path = os.path.join(temp_path,'NRGTEN', f'{target}_dynasig.pdb')
    with open(output_path, 'w') as t2:
        for line in texto_list:
            t2.write(line)
    return b_factor_dict


def prep_labels(labels):
    labels_list=[]
    for label in labels:
        res=int(label.split("|")[1])
        chain=label.split("|")[2]
        labels_list.append(f'{res}_{chain}')
    return(labels_list)


def dynamical_signature(form, install_dir):
    target = form.NRGten_select_target_object_1.currentText()
    lig = form.NRGten_select_ligand_object_1.currentText()
    target_2 = form.NRGten_select_target_object_2.currentText()
    beta = form.NRGten_dynasig_beta.text()
    main_folder_path = install_dir
    temp_path = form.temp_line_edit.text()

    start_time = time.time()
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    plots = []
    svib_list = []
    cmd.save(target_file, target)

    b_fact_dictionary_ref, dyna_sig_list_ref, model_ref, svib_ref = run_dynamical_signature(target_file, beta,
                                                                                            main_folder_path, temp_path)
    cmd.disable(target)
    selection = os.path.splitext(os.path.basename(target_file))[0] + '_dynasig'
    cmd.load(target_file[:-4] + '_dynasig.pdb')
    cmd.spectrum(selection=selection, palette='blue_white_red', expression='q')
    cmd.cartoon('putty', selection=selection)
    cmd.group('NRGTEN', selection)

    if target_2 == 'None':
        target_name = target
        svib_list.append(svib_ref)
        plots.append(go.Scatter(x=prep_labels(model_ref.get_mass_labels()), y=dyna_sig_list_ref, mode='lines',
                                name=f'Svib {svib_ref}'))

        if lig != 'None':
            output_file = os.path.join(temp_path,'NRGTEN', f'no_lig_{target}.pdb')
            remove_selection_and_save(target, lig, output_file)
            b_fact_dictionary_no_lig, dyna_sig_list_no_lig, model_no_lig, svib_no_lig = run_dynamical_signature(output_file, beta, main_folder_path, temp_path)
            svib_list.append(svib_no_lig)
            filename = os.path.splitext(os.path.basename(output_file))[0]

            for b_factor in range(len(dyna_sig_list_no_lig)):
                mass = model_no_lig.get_mass_labels()[b_factor]
                key = '{}_{}_{}'.format(mass.split('|')[0][:3], mass.split('|')[2], mass.split('|')[1])
                dyna_sig_list_no_lig[b_factor] = (b_fact_dictionary_ref[key] / dyna_sig_list_no_lig[b_factor]) - 1

            plots.append(go.Scatter(x=prep_labels(model_no_lig.get_mass_labels()), y=dyna_sig_list_no_lig, mode='lines',
                                    name=f'Svib {svib_no_lig}'))
            write_b_factor(filename, dyna_sig_list_no_lig, temp_path, model_no_lig.get_mass_labels())
            cmd.load(output_file[:-4] + '_dynasig.pdb')
            cmd.spectrum(selection=filename + '_dynasig', palette='blue_white_red', expression='q', minimum=-1,
                         maximum=1)
            cmd.cartoon('putty', selection=filename+ '_dynasig')
            cmd.group('NRGTEN', filename + '_dynasig')

    else:
        target_name = target_2
        object_list = []
        diff_list = []

        for state in range(cmd.count_states(target_2)):
            output_file = os.path.join(temp_path, 'NRGTEN', f'{target_2}_{state}.pdb')
            cmd.save(output_file, target_2, state=state + 1)

            diff = compare_residues(target_file, output_file)
            diff_list.append(diff)
            os.rename(output_file, os.path.join(temp_path, 'NRGTEN', f'{target_2}_{diff}.pdb'))
            output_file = os.path.join(temp_path, 'NRGTEN', f'{target_2}_{diff}.pdb')

            b_fact_dictionary_no_lig, dyna_sig_list_no_lig, model_no_lig, svib_no_lig = run_dynamical_signature(output_file, beta, main_folder_path, temp_path)
            svib_list.append(svib_no_lig - svib_ref)
            filename = os.path.splitext(os.path.basename(output_file))[0]

            for b_factor in range(len(dyna_sig_list_no_lig)):
                dyna_sig_list_no_lig[b_factor] = (dyna_sig_list_ref[b_factor] / dyna_sig_list_no_lig[b_factor]) - 1
            if 'LIG.' in model_no_lig.get_mass_labels()[-1]:
                dyna_sig_list_no_lig[-1] = 0
            dyna_sig_list_no_lig = standardize_to_minus1_plus1(dyna_sig_list_no_lig)

            plots.append(go.Scatter(x=prep_labels(model_no_lig.get_mass_labels()), y=dyna_sig_list_no_lig, mode='lines',
                                    name=f'Diff {diff}'), )
            write_b_factor(filename, dyna_sig_list_no_lig, temp_path, model_no_lig.get_mass_labels())
            cmd.load(os.path.join(temp_path, 'NRGTEN', f'{filename}_dynasig.pdb'), f'{target_2}_dynasigdif_{diff}')
            object_list.append(f'{target_2}_dynasigdif_{diff}')

        for state in diff_list:
            cmd.spectrum(selection=f'{target_2}_dynasigdif_{state}', palette='blue_white_red', expression='q',
                         minimum=-1, maximum=1)
            cmd.cartoon('putty', selection=f'{target_2}_dynasigdif_{state}')
        create_group(f'{target_2}_dynasigdif', object_list)
        cmd.group('NRGTEN', f'{target_2}_dynasigdif')

    fig = go.Figure()

    # Add traces but set them to be initially invisible, except for the first one
    for i, plot in enumerate(plots):
        fig.add_trace(plot)
        if i != 0:
            fig.data[i].visible = False

    # Create buttons to toggle visibility of each trace and for showing all plots together
    buttons = []
    all_visible_button = dict(label="All Combined", method="update",
                              args=[{"visible": [True for _ in range(len(plots))]}]
    )
    buttons.append(all_visible_button)

    for i in range(len(plots)):
        if target_2 == 'None':
            if i==1:
                label = f"Differential to target alone, DeltaSvib of target: {svib_list[i]:.2e}"
            else:
                label = f"Dynamical signature of complex, DeltaSvib: {svib_list[i]:.2e}"
        else:
            label = f"Diff {diff_list[i]}, DeltaSvib: {svib_list[i]:.2e}"
        button = dict(label=label, method="update", args=[{"visible": [j == i for j in range(len(plots))]}])
        buttons.append(button)

    # Update layout with the buttons
    fig.update_layout(
        updatemenus=[dict(type="buttons", showactive=True, buttons=buttons)],
        title=f"Dynamical Signatures of {target_name}",
        xaxis_title="Residue Index",
        yaxis_title="Fluctuation differential",
    )

    fig.write_html(os.path.join(temp_path, 'NRGTEN', f'{target_name}_diff.html'))
    fig.show()

    end_time = time.time()
    execution_time = end_time - start_time
    output_message(form.output_box, '=========== DynaSig ===========', 'valid')
    output_message(form.output_box, f"Execution time: {execution_time:.4f} seconds", 'valid')
    output_message(form.output_box, '=========== END DynaSig ===========', 'valid')

def create_group(group_name, object_list):
    members = ', '.join(object_list)
    cmd.group(group_name, members)
    return 0


def run_dynamical_signature(target_file, beta, main_folder_path, temp_path):
    _, filename = os.path.split(target_file)
    target = os.path.splitext(filename)[0]
    model = encom_model(target_file, main_folder_path, temp_path)
    dyna_sig = model.compute_bfactors_boltzmann(beta=float(beta))
    svib = model.compute_vib_entropy(beta=float(beta))
    b_fact_dictionary = write_b_factor(target, dyna_sig, temp_path, model.get_mass_labels())
    return b_fact_dictionary, dyna_sig, model, svib


def conformational_ensemble(form, install_dir):
    start_time = time.time()
    target = form.NRGten_select_target_object_1.currentText()
    modes_list = form.NRGten_modes_lineEdit.text()
    step = form.NRGten_step_lineEdit.text()
    max_conf = form.NRGten_max_conf_lineEdit.text()
    max_disp = form.NRGten_max_dis_lineEdit.text()
    opt = form.NRGten_optmizestates.isChecked()
    main_folder_path = install_dir
    temp_path = form.temp_line_edit.text()

    modes_list = list(map(int, modes_list.split(',')))
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    cmd.save(target_file, target)
    model = encom_model(target_file, main_folder_path, temp_path)
    ensemble_path = os.path.join(temp_path,'NRGTEN', f"{target}_conf_ensemble.pdb")
    model.build_conf_ensemble(modes_list, ensemble_path, step=float(step), max_displacement=float(max_disp),
                              max_conformations=int(max_conf))
    if opt:
        from srcs.nrgten.model_ensemble import model_states
        model_states(ensemble_path, target, temp_path, main_folder_path,form)

    else:
        cmd.load(ensemble_path)
        cmd.show('cartoon', f'{target}_conf_ensemble')
    end_time = time.time()
    execution_time = end_time - start_time
    output_message(form.output_box, '=========== Conformational Ensemble ===========', 'valid')
    output_message(form.output_box, f"Execution time: {execution_time:.4f} seconds", 'valid')
    output_message(form.output_box, '=========== END Conformational Ensemble ===========', 'valid')
