from general_functions import get_residue_info, process_flexaid_result
from PyQt5.QtCore import pyqtSignal, QThread
from modeller import Environ, Model, Selection, Alignment, log
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
from joblib import Parallel, delayed
import time
from general_functions import output_message
from pymol import cmd
import os
import multiprocessing
import general_functions
import subprocess
import sys
# TODO group residues by property

def flex_res(target_file):
    with open(target_file, "r") as f:
        texto=f.readlines()
        for line in texto:
            if 'LIG  9999' in line:
                return 1
    return 0

def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


#molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False



def check_all(form):
    checkbox_names = [
         'Modeller_checkBox', 'Modeller_checkBox_10', 'Modeller_checkBox_11',
         'Modeller_checkBox_12', 'Modeller_checkBox_13', 'Modeller_checkBox_14',
         'Modeller_checkBox_15', 'Modeller_checkBox_16', 'Modeller_checkBox_17',
         'Modeller_checkBox_18', 'Modeller_checkBox_19', 'Modeller_checkBox_2',
         'Modeller_checkBox_20', 'Modeller_checkBox_3', 'Modeller_checkBox_4',
         'Modeller_checkBox_5', 'Modeller_checkBox_6', 'Modeller_checkBox_7',
         'Modeller_checkBox_8', 'Modeller_checkBox_9'
    ]
    for checkbox_name in checkbox_names:
        checkbox = getattr(form, checkbox_name)
        if form.modeller_checkbox_all.isChecked():
            checkbox.setChecked(True)
        else:
            checkbox.setChecked(False)


class single_mutationManager:
    def __init__(self, form, install_dir):
        self.form = form
        self.target = form.Modeller_select_object.currentText()
        self.temp_path = form.temp_line_edit.text()
        self.selected_residue = self.form.Modeller_select_residue.currentText()
        self.install_dir = install_dir

    def run_single_mutation(self):
        cpu_usage_target = int(self.form.nrgrank_cpu_usage_target.currentText()[:-1])
        if cpu_usage_target == 100:
            number_of_cores = multiprocessing.cpu_count() + 4
        else:
            number_of_cores = round(multiprocessing.cpu_count() * (cpu_usage_target / 100))
        general_functions.disable_run_mutate_buttons(self.form, disable=True)
        residues_to_mutate = self.get_residues_to_mutate()
        residue_list = get_residue_info(self.selected_residue)
        number_of_mutations = len(residue_list) * len(residues_to_mutate)
        self.initialise_progress_bar(number_of_mutations)
        self.single_mutation_thread = single_mutationThread(self.temp_path, self.target, self.selected_residue, residues_to_mutate, number_of_cores, self.install_dir)
        self.single_mutation_thread.message_signal.connect(self.handle_message_signal)
        self.single_mutation_thread.mutation_progress_signal.connect(self.handle_mutation_progress_signal)
        self.single_mutation_thread.finished_signal.connect(self.handle_thread_finished)
        self.single_mutation_thread.start()

    def get_residues_to_mutate(self):
        checkbox_names = [
            'Modeller_checkBox', 'Modeller_checkBox_10', 'Modeller_checkBox_11',
            'Modeller_checkBox_12', 'Modeller_checkBox_13', 'Modeller_checkBox_14',
            'Modeller_checkBox_15', 'Modeller_checkBox_16', 'Modeller_checkBox_17',
            'Modeller_checkBox_18', 'Modeller_checkBox_19', 'Modeller_checkBox_2',
            'Modeller_checkBox_20', 'Modeller_checkBox_3', 'Modeller_checkBox_4',
            'Modeller_checkBox_5', 'Modeller_checkBox_6', 'Modeller_checkBox_7',
            'Modeller_checkBox_8', 'Modeller_checkBox_9'
        ]
        residues_to_mutate = []
        for checkbox_name in checkbox_names:
            checkbox = getattr(self.form, checkbox_name)
            if checkbox.isChecked():
                residues_to_mutate.append(checkbox.text())
        return residues_to_mutate

    def initialise_progress_bar(self, number_of_mutations):
        self.form.mutation_progress.show()
        self.form.mutation_progress.setEnabled(True)
        self.form.mutation_progress_label.setText(f'Mutation: 0/{number_of_mutations}')
        self.form.mutation_progress_bar.setValue(0)
        self.form.mutation_progress_bar.setMaximum(number_of_mutations)

    def handle_message_signal(self, message):
        general_functions.output_message(self.form.output_box, message, 'valid')

    def handle_mutation_progress_signal(self):
        progress_bar = self.form.mutation_progress_bar
        current_value = progress_bar.value()
        new_value = current_value + 1
        maximum_value = progress_bar.maximum()
        self.form.mutation_progress_bar.setValue(new_value)
        self.form.mutation_progress_label.setText(f'Mutation progress: {new_value}/{maximum_value}')

    def handle_thread_finished(self):
        self.single_mutation_thread.stop()
        self.single_mutation_thread.quit()
        self.single_mutation_thread.wait()
        self.single_mutation_thread = None
        general_functions.disable_run_mutate_buttons(self.form, enable=True)


class single_mutationThread(QThread):
    message_signal = pyqtSignal(str)
    mutation_progress_signal = pyqtSignal()
    finished_signal = pyqtSignal()

    def __init__(self, temp_path, target, selected_residue, desired_mutations, number_of_cores, install_dir):
        super().__init__()
        self.temp_path = temp_path
        self.target = target
        self.selected_residue = selected_residue
        self.desired_mutations = desired_mutations
        self.number_of_cores = number_of_cores
        self.install_dir = install_dir
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

    def run(self):
        start_time = time.time()
        self.message_signal.emit('=========== Single Mutations ===========')
        residue_selection = self.selected_residue
        target_file = os.path.join(self.temp_path, 'modeller', f'{self.target}.pdb')
        cmd.save(target_file, self.target)
        if flex_res(target_file):
            process_flexaid_result(target_file, target_file)
        residue_list = get_residue_info(residue_selection)
        mutant_files = []
        sorted_residue_list = sorted(residue_list, key=lambda x: (x[0], x[-1]))
        for residue in sorted_residue_list:
            residue_name = residue[0]
            residue_number = residue[1]
            residue_chain = residue[2]
            for desired_mutation in self.desired_mutations:
                self.mutation_progress_signal.emit()
                if desired_mutation != residue_name:
                    self.single_mutation_process = subprocess.Popen([sys.executable, os.path.join(self.install_dir, 'srcs', 'nrgsuite_modeller', 'modeller_one_mutation.py'), target_file, residue_number, desired_mutation, residue_chain])
                    while self.is_running and self.single_mutation_process.poll() is None:
                        self.msleep(100)
                    mutant_file = f'{target_file[:-4]}{desired_mutation}{residue_number}{residue_chain}.pdb'
                    mutant_files.append(mutant_file)
        for file_counter, file in enumerate(mutant_files):
            cmd.load(file, f'{self.target}_mutants', state=file_counter+1)
        cmd.group("Single_Mutants", self.target + '_mutants')
        cmd.disable(self.target)
        end_time = time.time()
        execution_time = end_time - start_time
        self.message_signal.emit(f"Execution time: {execution_time:.4f} seconds")
        self.message_signal.emit('=========== END Single Mutations ===========')
        self.finished_signal.emit()