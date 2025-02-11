import os
import shutil
from general_functions import disable_run_mutate_buttons, output_message
import subprocess
import sys
from PyQt5.QtCore import pyqtSignal, QThread

def delete_ligand_set(ligand_set_name, ligand_set_folder, output_box):
    if ligand_set_name == '':
        output_message(output_box, 'No ligand set name provided.', 'warning')
        return
    full_ligand_set_name = ligand_set_name.replace(' ', '_')
    ligand_set_path = os.path.join(ligand_set_folder, full_ligand_set_name)
    if os.path.exists(ligand_set_path):
        shutil.rmtree(ligand_set_path)
        output_message(output_box, f'Deleted ligand set {ligand_set_name}', 'valid')

class ConfGeneratorManager:
    def __init__(self, form, install_dir, ligand_set_folder_path):
        super().__init__()
        self.form = form
        self.install_dir = install_dir
        self.ligand_set_folder_path = ligand_set_folder_path
        self.new_smiles_file_path = None
        self.new_ligand_folder_path = None
        self.smiles_column_number = str(form.smiles_column_number.value())
        self.name_column_number = str(form.name_column_number.value())
        self.molecular_weight_max = str(form.molecular_weight_max.value())
        self.heavy_atoms_min = str(form.heavy_atoms_min.value())
        self.file_management()

    def file_management(self):
        smiles_path = self.form.nrgrank_add_ligand_file_path.text()
        new_ligand_folder_name = os.path.splitext(os.path.basename(smiles_path).replace(' ', '_'))[0]
        self.new_ligand_folder_path = os.path.join(self.ligand_set_folder_path, new_ligand_folder_name)
        if not os.path.exists(self.new_ligand_folder_path):
            os.makedirs(self.new_ligand_folder_path)
        self.new_smiles_file_path = os.path.join(self.new_ligand_folder_path, os.path.basename(smiles_path))
        shutil.copy(smiles_path, self.new_smiles_file_path)

    def generate_conformer(self):
        self.generate_conformer_thread = GenerateConformerThread(self.install_dir,
                                                                 self.new_smiles_file_path,
                                                                 self.new_ligand_folder_path,
                                                                 self.smiles_column_number,
                                                                 self.name_column_number,
                                                                 self.molecular_weight_max,
                                                                 self.heavy_atoms_min)
        self.generate_conformer_thread.message_signal.connect(self.handle_message_signal)
        self.generate_conformer_thread.finished_signal.connect(self.handle_thread_finished)
        self.generate_conformer_thread.start()

    def handle_message_signal(self, message, message_type):
        output_message(self.form.output_box, message, message_type)

    def handle_thread_finished(self):
        self.generate_conformer_thread.stop()
        self.generate_conformer_thread.quit()
        self.generate_conformer_thread.wait()
        self.nrgrank_thread = None
        disable_run_mutate_buttons(self.form, enable=True)


class GenerateConformerThread(QThread):
    message_signal = pyqtSignal(str,str)
    finished_signal = pyqtSignal()

    def __init__(self, install_dir, smiles_path, custom_output_path, smiles_column_number, name_column_number,
                 molecular_weight_max, heavy_atoms_min):
        super().__init__()
        self.install_dir = install_dir
        self.smiles_path = smiles_path
        self.ligand_set_path = custom_output_path
        self.is_running = True
        self.deps_path = os.path.join(self.install_dir, 'deps', 'nrgrank')
        self.config_path = os.path.join(self.deps_path, 'config.txt')
        self.smiles_column_number = smiles_column_number
        self.name_column_number = name_column_number
        self.molecular_weight_max = molecular_weight_max
        self.heavy_atoms_min = heavy_atoms_min

    def stop(self):
        self.is_running = False

    def run(self):
        try:
            self.message_signal.emit('Starting to generate conformers', 'valid')
            generate_conformer_path = os.path.join(self.install_dir, 'srcs', 'nrgrank', 'generate_conformers.py')
            conformer_folder_path = os.path.splitext(self.smiles_path)[0] + '_conformers'
            smile_file_name = os.path.splitext(os.path.basename(self.smiles_path))[0]
            mol2_path = os.path.join(conformer_folder_path, smile_file_name + '_1_conf.mol2')
            preprocessed_ligand_path = os.path.join(conformer_folder_path, 'preprocessed_ligands_1_conf')

            self.generate_conformer_process = subprocess.Popen(
                [sys.executable, generate_conformer_path, '-s', self.smiles_path, '-d', self.deps_path,
                 '-sc', self.smiles_column_number, '-nc', self.name_column_number, '-mw', self.molecular_weight_max, '-ha', self.heavy_atoms_min],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            while self.is_running and self.generate_conformer_process.poll() is None:
                self.msleep(100)


            stdout, stderr = self.generate_conformer_process.communicate()

            if self.generate_conformer_process.returncode != 0:
                self.message_signal.emit(f"Error while generating conformer:\n{stderr}", 'error')
            if self.generate_conformer_process.returncode == 0:
                output_file_path = os.path.join(self.ligand_set_path, 'conformer_generation.out')
                with open(output_file_path, 'w') as output_file:
                    output_file.write(stdout)
            if os.path.isdir(preprocessed_ligand_path):
                shutil.move(preprocessed_ligand_path, self.ligand_set_path)
            if os.path.isfile(mol2_path):
                shutil.move(mol2_path, self.ligand_set_path)
            if os.path.isdir(conformer_folder_path):
                shutil.rmtree(conformer_folder_path)
            self.message_signal.emit(f"Finished generating conformers", 'valid')
            self.finished_signal.emit()


        except Exception as e:
            self.message_signal.emit(f"An unexpected error occurred: {e}", 'error')
