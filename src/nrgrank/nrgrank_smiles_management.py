import os
import shutil
from general_functions import output_message
import subprocess
import sys
from PyQt5.QtCore import pyqtSignal, QThread

def delete_ligand_set(ligand_set_name, ligand_set_folder, output_box):
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
        print('here')
        self.generate_conformer_thread = GenerateConformerThread(self.install_dir, self.new_smiles_file_path, 'no', self.new_ligand_folder_path)
        self.generate_conformer_thread.start()

class GenerateConformerThread(QThread):

    def __init__(self, install_dir, smiles_path, optimize, custom_output_path):
        super().__init__()
        self.install_dir = install_dir
        self.smiles_path = smiles_path
        self.optimize = optimize
        self.custom_output_path = custom_output_path
        self.is_running = True
        self.deps_path = os.path.join(self.install_dir, 'deps', 'nrgrank')
        self.config_path = os.path.join(self.deps_path, 'config.txt')

    def run(self):
        print('running NRGRank thread')
        generate_conformer_path = os.path.join(self.install_dir, 'src', 'nrgrank','generate_conformers.py')
        sdf_path = os.path.splitext(self.smiles_path)[0] + '.sdf'
        mol2_path = os.path.splitext(self.smiles_path)[0] + '.mol2'
        self.generate_conformer_process = subprocess.Popen([sys.executable,
                                                        generate_conformer_path,
                                                        self.smiles_path,
                                                        self.optimize,
                                                        sdf_path,
                                                        self.deps_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while self.is_running and self.generate_conformer_process.poll() is None:
            self.msleep(100)
        if os.path.isfile(mol2_path):
            os.remove(mol2_path)