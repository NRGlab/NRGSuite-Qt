from general_functions import get_residue_info, process_flexaid_result
from modeller import Environ, Model, Selection, Alignment, log
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import pyqtSignal, QThread
from joblib import Parallel, delayed
import time
from general_functions import output_message
from pymol import cmd
import os
import multiprocessing
import general_functions
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


#use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = Selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)


def sm(modelname, respos, restyp, chain, temp_path):
    log.none()

    # Set a different value for rand_seed to get a different final model
    env = Environ(rand_seed=-49837)

    env.io.hetatm = True
    #soft sphere potential
    env.edat.dynamic_sphere=False
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = Model(env, file=modelname)

    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    #set up the mutate residue selection segment
    s = Selection(mdl1.chains[chain].residues[respos])

    #perform the mutate residue operation
    s.mutate(residue_type=restyp)
    #get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])


    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    #here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    #yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = Model(env, file=modelname)

    #required to do a transfer_res_numb
    #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    #transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

    #It is usually necessary to write the mutated sequence out and read it in
    #before proceeding, because not all sequence related information about MODEL
    #is changed by this command (e.g., internal coordinates, charges, and atom
    #types and radii are not updated).

    mdl1.write(file=modelname+restyp+respos+'.tmp')
    mdl1.read(file=modelname+restyp+respos+'.tmp')

    #set up restraints before computing energy
    #we do this a second time because the model has been written out and read in,
    #clearing the previously set restraints
    make_restraints(mdl1, ali)

    #a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms=1

    sched = autosched.loop.make_for_model(mdl1)

    #only optimize the selected residue (in first pass, just atoms in selected
    #residue, in second pass, include nonbonded neighboring atoms)
    #set up the mutate residue selection segment
    s = Selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()

    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

    #feels environment (energy computed on pairs that have at least one member
    #in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)

    s.energy()
    #give a proper name
    mdl1.write(modelname[:-4]+restyp+respos+'.pdb')
    #delete the temporary file
    os.remove(modelname+restyp+respos+'.tmp')
    return 0


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


def model_mutations(form, temp_path):
    start_time = time.time()
    target = form.Modeller_select_object.currentText()
    res_list = form.Modeller_select_residue.currentText()
    mutation_list = [form.Modeller_checkBox.isChecked(),
                     form.Modeller_checkBox_2.isChecked(),
                     form.Modeller_checkBox_3.isChecked(),
                     form.Modeller_checkBox_4.isChecked(),
                     form.Modeller_checkBox_5.isChecked(),
                     form.Modeller_checkBox_6.isChecked(),
                     form.Modeller_checkBox_7.isChecked(),
                     form.Modeller_checkBox_8.isChecked(),
                     form.Modeller_checkBox_9.isChecked(),
                     form.Modeller_checkBox_10.isChecked(),
                     form.Modeller_checkBox_11.isChecked(),
                     form.Modeller_checkBox_12.isChecked(),
                     form.Modeller_checkBox_13.isChecked(),
                     form.Modeller_checkBox_14.isChecked(),
                     form.Modeller_checkBox_15.isChecked(),
                     form.Modeller_checkBox_16.isChecked(),
                     form.Modeller_checkBox_17.isChecked(),
                     form.Modeller_checkBox_18.isChecked(),
                     form.Modeller_checkBox_19.isChecked(),
                     form.Modeller_checkBox_20.isChecked()]
    target_file = os.path.join(temp_path,'modeller', f'{target}.pdb')
    cmd.save(target_file, target)
    if flex_res(target_file):
        process_flexaid_result(target_file, target_file)
    res_list=get_residue_info(res_list)
    amino_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    count=1
    for res in res_list:
        for res_1 in range(len(amino_list)):
            if mutation_list[res_1]:
                if amino_list[res_1] != res[0]:
                        sm(target_file,res[1],amino_list[res_1],res[2],temp_path)
                        cmd.load(target_file[:-4]+amino_list[res_1]+res[1]+'.pdb', target + '_mutants',state=count)
                        count+=1

    cmd.group("Single_Mutants", target + '_mutants')
    end_time = time.time()
    execution_time = end_time - start_time
    output_message(form.output_box, '=========== Single Mutations ===========', 'valid')
    output_message(form.output_box, f"Execution time: {execution_time:.4f} seconds", 'valid')
    output_message(form.output_box, '=========== END Single Mutations ===========', 'valid')


class single_mutationManager:
    def __init__(self, form):
        self.form = form
        self.target = form.Modeller_select_object.currentText()
        self.temp_path = form.temp_line_edit.text()
        self.selected_residue = self.form.Modeller_select_residue.currentText()

    def run_single_mutation(self):
        cpu_usage_target = int(self.form.nrgrank_cpu_usage_target.currentText()[:-1])
        if cpu_usage_target == 100:
            number_of_cores = multiprocessing.cpu_count() + 4
        else:
            number_of_cores = round(multiprocessing.cpu_count() * (cpu_usage_target / 100))
        general_functions.disable_run_mutate_buttons(self.form, disable=True)
        residues_to_mutate = self.get_residues_to_mutate()
        self.single_mutation_thread = single_mutationThread(self.temp_path, self.target, self.selected_residue, residues_to_mutate, number_of_cores)
        self.single_mutation_thread.message_signal.connect(self.handle_message_signal)
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

    def handle_message_signal(self, message):
        general_functions.output_message(self.form.output_box, message, 'valid')

    def handle_screen_progress_signal(self, value):
        current_value = self.form.nrgrank_progress_bar.value()
        if value > current_value:
            self.form.nrgrank_progress_bar.setValue(value)
            self.form.nrgrank_progress_label.setText(f'Screening progress: {value}%')

    def handle_thread_finished(self):
        self.single_mutation_thread.stop()
        self.single_mutation_thread.quit()
        self.single_mutation_thread.wait()
        self.single_mutation_thread = None
        general_functions.disable_run_mutate_buttons(self.form, enable=True)


class single_mutationThread(QThread):
    message_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, temp_path, target, selected_residue, desired_mutations, number_of_cores):
        super().__init__()
        self.temp_path = temp_path
        self.target = target
        self.selected_residue = selected_residue
        self.desired_mutations = desired_mutations
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
        for residue in residue_list:
            residue_name = residue[0]
            residue_number = residue[1]
            residue_chain = residue[2]
            for desired_mutation in self.desired_mutations:
                if desired_mutation != residue_name:
                    print('Desired mutation:', desired_mutation)
                    sm(target_file, residue_number, desired_mutation, residue_chain, self.temp_path)
                    mutant_file = f'{target_file[:-4]}{desired_mutation}{residue_number}.pdb'
                    mutant_files.append(mutant_file)
        for file_counter, file in enumerate(mutant_files):
            cmd.load(file, f'{self.target}_mutants', state=file_counter+1)

        cmd.group("Single_Mutants", self.target + '_mutants')
        end_time = time.time()
        execution_time = end_time - start_time
        self.message_signal.emit(f"Execution time: {execution_time:.4f} seconds")
        self.message_signal.emit('=========== END Single Mutations ===========')
        self.finished_signal.emit()