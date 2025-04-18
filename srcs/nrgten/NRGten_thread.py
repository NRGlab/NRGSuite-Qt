import general_functions
from PyQt5.QtCore import pyqtSignal, QThread
import time, os, multiprocessing, json, sys, pickle
import numpy as np
from pymol import cmd
import plotly.graph_objects as go
from srcs.surfaces.ligand_atomtypes import add_pdb
from srcs.surfaces.run_Surfaces import create_ligand_file, compare_residues
from nrgten.atom import Atom
import subprocess
from joblib import Parallel, delayed
from general_functions import process_flexaid_result


class DynasigManager:
    def __init__(self, form, install_dir):
        self.form = form
        self.target = form.NRGten_select_target_object_1.currentText()
        self.lig = form.NRGten_select_ligand_object_1.currentText()
        self.target_2 = form.NRGten_select_target_object_2.currentText()
        self.beta = form.NRGten_dynasig_beta.text()
        self.main_folder_path = install_dir
        self.temp_path = form.temp_line_edit.text()

    def run_nrgten(self):
        cpu_usage_target = int(self.form.nrgrank_cpu_usage_target.currentText()[:-1])
        if cpu_usage_target == 100:
            number_of_cores = multiprocessing.cpu_count() + 4
        else:
            number_of_cores = round(multiprocessing.cpu_count() * (cpu_usage_target / 100))

        general_functions.disable_run_mutate_buttons(self.form, disable=True)
        self.dynasig_thread = DynasigThread(self.temp_path, self.target, self.beta, self.lig, self.target_2, self.main_folder_path, number_of_cores)
        self.dynasig_thread.message_signal.connect(self.handle_message_signal)
        self.dynasig_thread.finished_signal.connect(self.handle_thread_finished)
        self.dynasig_thread.start()

    def handle_message_signal(self, message):
        general_functions.output_message(self.form.output_box, message, 'valid')

    def handle_screen_progress_signal(self, value):
        current_value = self.form.nrgrank_progress_bar.value()
        if value > current_value:
            self.form.nrgrank_progress_bar.setValue(value)
            self.form.nrgrank_progress_label.setText(f'Screening progress: {value}%')

    def handle_thread_finished(self):
        self.dynasig_thread.stop()
        self.dynasig_thread.quit()
        self.dynasig_thread.wait()
        self.dynasig_thread = None
        general_functions.disable_run_mutate_buttons(self.form, enable=True)


class DynasigThread(QThread):
    message_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, temp_path, target, beta, lig, target_2, main_folder_path, number_of_cores):
        super().__init__()
        self.temp_path = temp_path
        self.nrgten_temp_path = os.path.join(temp_path, 'NRGTEN')
        self.target = target
        self.beta = beta
        self.lig = lig
        self.target_2 = target_2
        self.main_folder_path = main_folder_path
        self.number_of_cores = number_of_cores
        self.is_running = True
        self.processes = []
        self.executor = None

    def stop(self):
        self.is_running = False
        if self.executor:
            self.executor.shutdown(wait=False)
        for process in self.processes:
            if process.poll() is None:
                process.terminate()
                process.wait()

    @staticmethod
    def prep_labels(labels):
        labels_list = []
        for label in labels:
            res = int(label.split("|")[1])
            chain = label.split("|")[2]
            labels_list.append(f'{res}_{chain}')
        return (labels_list)

    @staticmethod
    def remove_selection_and_save(object_name, selection, output_file):
        cmd.create("object_without_selection", f"{object_name} and not ({selection})")
        cmd.save(output_file, "object_without_selection")
        cmd.delete("object_without_selection")


    @staticmethod
    def standardize_to_minus1_plus1(data):
        max_abs_value = max(abs(x) for x in data)
        standardized_data = [x / max_abs_value for x in data]
        return standardized_data

    @staticmethod
    def create_group(group_name, object_list):
        members = ', '.join(object_list)
        cmd.group(group_name, members)
        return 0

    @staticmethod
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

    def find_het(self, target_file, main_folder_path):
        het_dic = {}
        with open(target_file, "r") as t1:
            texto = t1.readlines()
            for line in texto:
                if 'HETATM' in line:
                    het_dic[line[17:20]] = '1'
        for lig in list(het_dic):
            def_file = os.path.join(main_folder_path, "deps", "surfaces", 'AMINO_FlexAID.def')
            open_def_file = open(def_file, "r")
            ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
            create_ligand_file(target_file, ligand_file_name)
            custom_def_path = os.path.join(self.nrgten_temp_path, f'custom_{os.path.basename(def_file)}')
            custom_def_file = open(custom_def_path, 'w')
            add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
            open_def_file.close()
            custom_def_file.close()
            self.generate_massfile(ligand_file_name + '.pdb', os.path.join(self.nrgten_temp_path, f'{lig}.masses'))
            with open(custom_def_path, "r") as t2:
                texto = t2.readlines()
                definition = texto[-1][:-2] + '\n'
            with open(os.path.join(self.nrgten_temp_path, f'{lig}.atomtypes'), "w") as t3:
                t3.write(definition)
        return list(het_dic)

    @staticmethod
    def write_b_factor(target, dyna_sig, nrgten_temp_path, labels):
        target_file = os.path.join(nrgten_temp_path, f'{target}.pdb')
        b_factor_dict = {}
        with open(os.path.join(nrgten_temp_path, target + '_dynasig.txt'), 'w') as t1:
            for res in range(len(labels)):
                t1.write('{} {}\n'.format(labels[res], dyna_sig[res]))
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
                    lined_abs = f"{abs(number):.2f}"
                    lined = lined.rjust(6)
                    lined_abs = lined_abs.rjust(6)
                    texto_list.append(line[:54] + lined + lined_abs + line[66:])
                else:
                    texto_list.append(line)
        output_path = os.path.join(nrgten_temp_path, f'{target}_dynasig.pdb')
        with open(output_path, 'w') as t2:
            for line in texto_list:
                t2.write(line)
        return b_factor_dict

    @staticmethod
    def process_state(state, state_pdb_file, list_het, temp_path, main_folder_path, beta):
        print(f'= State {state} started =')
        command = [sys.executable,
                   os.path.join(main_folder_path, 'srcs', 'nrgten', 'nrgten_separate.py'),
                   '-t', state_pdb_file,
                   '-b', beta,
                   '-m', main_folder_path,
                   '-te', temp_path,
                   '-l', list_het]
        with subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                              close_fds=True) as proc:
            proc.wait()
        pickle_file_path = os.path.splitext(state_pdb_file)[0] + '.pkl'
        with open(pickle_file_path, "rb") as f:
            b_fact_dictionary_no_lig, dyna_sig_list_no_lig, model_no_lig_mass_label, svib_no_lig = pickle.load(f)
        return b_fact_dictionary_no_lig, dyna_sig_list_no_lig, model_no_lig_mass_label, svib_no_lig

    def run(self):
        start_time = time.time()
        self.message_signal.emit("=========== DynaSig ===========")
        target_file = os.path.join(self.nrgten_temp_path, f'{self.target}.pdb')
        plots = []
        svib_list = []
        cmd.save(target_file, self.target)
        process_flexaid_result(target_file,target_file)
        list_het = self.find_het(target_file, self.main_folder_path)
        list_het = json.dumps(list_het)
        process_init_dynamic_signature = subprocess.Popen([sys.executable,
                                                        os.path.join(self.main_folder_path, 'srcs', 'nrgten',
                                                                     'nrgten_separate.py'),
                                                        '-t', target_file,
                                                        '-b', self.beta,
                                                        '-m', self.main_folder_path,
                                                        '-te', self.temp_path,
                                                        '-l', list_het])
        while self.is_running and process_init_dynamic_signature.poll() is None:
            self.msleep(100)
        pickle_file_path = os.path.splitext(target_file)[0] + '.pkl'
        with open(pickle_file_path, "rb") as f:
            b_fact_dictionary_ref, dyna_sig_list_ref, model_ref_mass_labels, svib_ref = pickle.load(f)
        cmd.disable(self.target)
        selection = os.path.splitext(os.path.basename(target_file))[0] + '_dynasig'
        cmd.load(target_file[:-4] + '_dynasig.pdb')
        cmd.spectrum(selection=selection, palette='blue_white_red', expression='q')
        cmd.cartoon('putty', selection=selection)
        cmd.group('NRGTEN', selection)
        if self.target_2 == 'None':
            target_name = self.target
            svib_list.append(svib_ref)
            plots.append(go.Scatter(x=self.prep_labels(model_ref_mass_labels), y=dyna_sig_list_ref, mode='lines',
                                    name=f'Svib {svib_ref}'))

            if self.lig != 'None':
                output_file = os.path.join(self.nrgten_temp_path, f'no_lig_{self.target}.pdb')
                self.remove_selection_and_save(self.target, self.lig, output_file)
                list_het_no_lig = self.find_het(output_file, self.main_folder_path)
                list_het_no_lig = json.dumps(list_het_no_lig)
                p = subprocess.Popen([sys.executable, os.path.join(self.main_folder_path, 'srcs', 'nrgten','nrgten_separate.py'),
                                                                   '-t', output_file,
                                                                   '-b', self.beta,
                                                                   '-m', self.main_folder_path,
                                                                   '-te', self.temp_path,
                                                                   '-l', list_het_no_lig])
                while self.is_running and p.poll() is None:
                    self.msleep(100)
                pickle_file_path = os.path.splitext(output_file)[0] + '.pkl'
                with open(pickle_file_path, "rb") as f:
                    b_fact_dictionary_no_lig, dyna_sig_list_no_lig, model_no_lig_mass_label, svib_no_lig = pickle.load(f)
                svib_list.append(svib_no_lig)
                filename = os.path.splitext(os.path.basename(output_file))[0]

                for b_factor in range(len(dyna_sig_list_no_lig)):
                    mass = model_no_lig_mass_label[b_factor]
                    key = '{}_{}_{}'.format(mass.split('|')[0][:3], mass.split('|')[2], mass.split('|')[1])
                    dyna_sig_list_no_lig[b_factor] = (b_fact_dictionary_ref[key] / dyna_sig_list_no_lig[b_factor]) - 1

                plots.append(
                    go.Scatter(x=self.prep_labels(model_no_lig_mass_label), y=dyna_sig_list_no_lig, mode='lines',
                               name=f'Svib {svib_no_lig}'))
                self.write_b_factor(filename, dyna_sig_list_no_lig, self.nrgten_temp_path, model_no_lig_mass_label)
                cmd.load(output_file[:-4] + '_dynasig.pdb')
                cmd.spectrum(selection=filename + '_dynasig', palette='blue_white_red', expression='q', minimum=-1,
                             maximum=1)
                cmd.cartoon('putty', selection=filename + '_dynasig')
                cmd.group('NRGTEN', filename + '_dynasig')

        else:
            target_name = self.target_2
            self.message_signal.emit(f"Starting states")
            states = range(cmd.count_states(self.target_2))
            state_file_list = []
            diff_list = []
            object_list = []
            for state_counter, state in enumerate(states):
                output_file_state = os.path.join(self.nrgten_temp_path, f'{self.target_2}_{state}.pdb')
                cmd.save(output_file_state, self.target_2, state=state + 1)

                diff = compare_residues(target_file, output_file_state)
                output_file_diff = os.path.join(self.nrgten_temp_path, f'{self.target_2}_{diff}.pdb')
                if os.path.isfile(output_file_diff):
                    os.remove(output_file_diff)
                os.rename(output_file_state, output_file_diff)
                state_file_list.append(output_file_diff)
                diff_list.append(diff)
                object_name = f'{self.target_2}_dynasigdif_{diff_list[state_counter]}'
                object_list.append(object_name)
            results = Parallel(n_jobs=-1)(delayed(self.process_state)(state, state_file_list[state],
                                                                     list_het, self.temp_path, self.main_folder_path,
                                                                     self.beta) for state in states)
            b_fact_dictionary_list_no_lig = [result[0] for result in results]
            dyna_sig_list_list_no_lig = [result[1] for result in results]
            model_list_no_lig_mass_label = [result[2] for result in results]
            svib_list_no_lig = [result[3] for result in results]

            for state in states:
                state_file= state_file_list[state]
                _ = b_fact_dictionary_list_no_lig[state]
                dyna_sig_list_no_lig = dyna_sig_list_list_no_lig[state]
                model_no_lig_mass_label = model_list_no_lig_mass_label[state]
                svib_no_lig = svib_list_no_lig[state]

                svib_list.append(svib_no_lig - svib_ref)
                for b_factor in range(len(dyna_sig_list_no_lig)):
                    dyna_sig_list_no_lig[b_factor] = (dyna_sig_list_ref[b_factor] / dyna_sig_list_no_lig[b_factor]) - 1
                if 'LIG.' in model_no_lig_mass_label[-1]:
                    dyna_sig_list_no_lig[-1] = 0

                dyna_sig_list_no_lig = self.standardize_to_minus1_plus1(dyna_sig_list_no_lig)

                filename = os.path.splitext(os.path.basename(state_file))[0]
                plot = go.Scatter(x=self.prep_labels(model_no_lig_mass_label),
                                  y=dyna_sig_list_no_lig, mode='lines',
                                  name=f'Diff {diff_list[state]}')
                plots.append(plot)
                self.write_b_factor(filename, dyna_sig_list_no_lig, self.nrgten_temp_path, model_no_lig_mass_label)
                cmd.load(os.path.join(self.nrgten_temp_path, f'{filename}_dynasig.pdb'), object_list[state])
                cmd.spectrum(selection=object_list[state], palette='blue_white_red', expression='q',
                             minimum=-1, maximum=1)
                cmd.cartoon('putty', selection=object_list[state])
            self.create_group(f'{self.target_2}_dynasigdif', object_list)
            cmd.group('NRGTEN', f'{self.target_2}_dynasigdif')
        if self.is_running:
            fig = go.Figure()

            for i, plot in enumerate(plots):
                fig.add_trace(plot)
                if i != 0:
                    fig.data[i].visible = False

            buttons = []
            all_visible_button = dict(label="All Combined", method="update",
                                      args=[{"visible": [True for _ in range(len(plots))]}]
                                      )
            buttons.append(all_visible_button)

            for i in range(len(plots)):
                if self.target_2 == 'None':
                    if i == 1:
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

            fig.write_html(os.path.join(self.nrgten_temp_path, f'{target_name}_diff.html'))
            fig.show()

            end_time = time.time()
            execution_time = end_time - start_time
            self.message_signal.emit(f"Execution time: {execution_time:.4f} seconds")
            self.message_signal.emit('=========== END DynaSig ===========')
            self.finished_signal.emit()

