from general_functions import get_residue_info, process_flexaid_result
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
from PyQt5.QtWidgets import QWidget


from pymol import cmd
import os
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
    log.verbose()

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
    objects=form.findChildren(QWidget)
    checkbox_names = [
         'Modeller_checkBox', 'Modeller_checkBox_10', 'Modeller_checkBox_11',
         'Modeller_checkBox_12', 'Modeller_checkBox_13', 'Modeller_checkBox_14',
         'Modeller_checkBox_15', 'Modeller_checkBox_16', 'Modeller_checkBox_17',
         'Modeller_checkBox_18', 'Modeller_checkBox_19', 'Modeller_checkBox_2',
         'Modeller_checkBox_20', 'Modeller_checkBox_3', 'Modeller_checkBox_4',
         'Modeller_checkBox_5', 'Modeller_checkBox_6', 'Modeller_checkBox_7',
         'Modeller_checkBox_8', 'Modeller_checkBox_9'
    ]
    for obj in objects:
        if obj.objectName() in checkbox_names:
            obj.setChecked(not(obj.isChecked()))


# Example Usage
def model_mutations(form, temp_path):
    # Specify input PDB file, residue number, and new residue
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
    target_file = os.path.join(temp_path,'modeller', '{}.pdb'.format(target))
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
                        cmd.load(target_file[:-4]+amino_list[res_1]+res[1]+'.pdb',target+'_mutants',state=count)
                        count+=1
    cmd.group("Single Mutants",target+'_mutants')

