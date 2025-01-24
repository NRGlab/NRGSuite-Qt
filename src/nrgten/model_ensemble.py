from modeller import Environ, Model, Alignment
from modeller.automodel import *
from modeller.selection import selection
from pymol import cmd
from src.surfaces.clean_structure import main as clean_structure
from PyQt5.QtWidgets import QApplication
import os

def model_states(input_file, target, temp_path, main_folder_path, form):
    model_data = []
    current_model = []
    model_dir = os.path.join(temp_path,'NRGTEN', f'{target}_models')
    os.makedirs(model_dir,exist_ok=True)

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('MODEL'):
                if current_model:
                    model_data.append(current_model)
                current_model = [line]
            elif line.startswith('ENDMDL'):
                current_model.append(line)
                model_data.append(current_model)
                current_model = []
            else:
                if current_model:
                    current_model.append(line)
    if current_model:
        model_data.append(current_model)

    # Set the maximum value for the progress bar
    form.NRGten_progress.setMaximum(len(model_data))

    # Iterate over each model and process it
    for i, model in enumerate(model_data):
        model_output_path = os.path.join(model_dir, f'model_{i + 1}.pdb')
        with open(model_output_path, 'w') as output_file:
            output_file.writelines(model)

        code = os.path.join(model_dir, f'model_{i + 1}')
        cleaned_file_path = code + '_clean.pdb'

        # Clean PDB file
        clean_pdb_file = open(cleaned_file_path, "w")
        clean_structure(code + '.pdb', os.path.join(main_folder_path, "deps", "surfaces", 'AMINO_FlexAID.def'),
                        clean_pdb_file)
        clean_pdb_file.close()

        # Run model state alignment and generation
        model_state(cleaned_file_path, model_dir)

        # Update UI elements: progress bar and label
        form.NRGten_progress.setValue(i + 1)
        form.NRGten_label.setText(f'State: {i + 1}/{len(model_data)}')

        # Process events to keep the UI responsive
        QApplication.processEvents()

        # Load model into PyMOL and show as cartoon
        cmd.load(cleaned_file_path, f'{target}_ensemble', state=i + 1)
        cmd.show('cartoon', f'{target}_ensemble')



def model_state(cleaned_file_path, model_dir):
    pdb_path = cleaned_file_path
    env = Environ()
    mdl = Model(env, file=pdb_path)
    aln = Alignment(env)
    aln.append_model(mdl, align_codes='pdb_target', atom_files=pdb_path)
    mdl = Model(env, file=pdb_path)
    aln.append_model(mdl, align_codes='pdb_template', atom_files=pdb_path)
    aln.align2d()
    aln.write(file='alignment.ali', alignment_format='PIR')

    class MyModel(automodel):
        def select_atoms(self):
            return selection(self)

    output_dir = os.path.join(model_dir, 'output_models')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    a = MyModel(env, alnfile='alignment.ali', knowns='pdb_template', sequence='pdb_target')
    a.starting_model = 1
    a.ending_model = 1
    a.very_fast()
    a.make()
    os.rename('pdb_target.B99990001.pdb', cleaned_file_path)