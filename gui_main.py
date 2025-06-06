import os
import sys

install_dir = os.path.dirname(__file__)
sys.path.append(install_dir)

import shutil
import subprocess
import platform
import general_functions
from srcs.flexaid.flexaid import FlexAIDManager, stop_simulation, abort_simulation, pause_resume_simulation, flexaid_show_ligand_from_table
from srcs.getcleft import getcleft
from srcs.nrgrank import nrgrank_on_target
from srcs.getcleft import spheres
from srcs.surfaces import run_Surfaces
from srcs.isomif import run_isomif
from srcs.nrgrank import nrgrank_smiles_management
from srcs.nrgten.NRGten_thread import DynasigManager
from srcs.settings import run_settings
from PyQt5.QtWidgets import QWidget
from PyQt5.QtGui import QStandardItemModel
from PyQt5.uic import loadUi
from srcs.nrgten import run_NRGTEN
try:
    import modeller
except ImportError:
    print('Modeller not installed.')
else:
    from srcs.nrgsuite_modeller import run_modeller
    from srcs.nrgsuite_modeller.run_modeller import single_mutationManager
# TODO: when showing surfaces result hide everything else
# TODO: clickable results in nrgrank table

def test_binary(binary_folder_path, operating_system):
    all_files = os.listdir(binary_folder_path)
    binary_files = []
    for f in all_files:
        full_path = os.path.join(binary_folder_path, f)
        if os.path.isfile(full_path) and not f.startswith('.'):
            binary_files.append(full_path)
    if operating_system == 'mac':
        for file in binary_files:
            subprocess.run(["chmod", "755", file])
            if file.endswith('isomif'):
                result = subprocess.run([file, '-h'], capture_output=True, text=True)
            else:
                result = subprocess.run([file], capture_output=True, text=True)
            if result.returncode != 0 and result.returncode != -11 and result.returncode != 24:
                print('Could not run: ', file)


class Controller:
    def __init__(self, form, binary_folder_path, binary_suffix, operating_system, ligand_set_folder_path, color_list):
        self.form = form
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.operating_system = operating_system
        self.ligand_set_folder_path = ligand_set_folder_path
        self.color_list = color_list
        # self.form.nrgrank_result_table.setSelectionMode(QTableWidget.MultiSelection)
        self.model = QStandardItemModel()
        self.form.nrgrank_result_table.setModel(self.model)
        self.form.flexaid_result_table.setModel(self.model)
        self.setupConnections()
        self.form.stackedWidget.setCurrentIndex(7)
        self.advanced_settings_dialog = None

    def setupConnections(self):
        self.form.button_getcleft.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(0))
        self.form.button_nrgrank.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(1))
        self.form.button_flexaid.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(2))
        self.form.button_surfaces.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(3))
        self.form.button_nrgten.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(4))
        self.form.button_modeller.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(5))
        self.form.button_ISOMIF.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(6))
        self.form.button_home.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(7))
        self.form.button_settings.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(8))
        self.form.button_tools.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(9))


        # Save/load session
        self.form.button_save.clicked.connect(lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text()))
        self.form.button_load.clicked.connect(lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text(), save=0))

        # GetCleft
        self.form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(self.form))
        self.form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box))
        self.form.button_start.clicked.connect(self.run_getcleft)
        #TODO: colour with cleft and name to be able to delete not important ones

        # Partition Cleft
        self.form.cleft_partition_button_add.clicked.connect(
            lambda: spheres.display_sphere(self.form.cleft_partition_select_object.currentText(), self.form,
                                           self.form.cleft_partition_radius_slider, self.form.temp_line_edit.text()))
        self.form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.cleft_partition_select_object, self.form.output_box,filter_for='bd_site'))
        self.form.cleft_partition_button_move.clicked.connect(spheres.move_sphere)
        self.form.cleft_partition_radius_slider.valueChanged.connect( lambda: spheres.resize_sphere('SPHERE', self.form.cleft_partition_radius_slider.value()))
        self.form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft('SPHERE', self.form.cleft_partition_radius_slider.value() / 100, self.form.temp_line_edit.text(), self.form.cleft_partition_select_object.currentText(), self.form.output_box, self.form.cleft_partition_radius_slider))
        self.form.cleft_partition_button_delete.clicked.connect(lambda: spheres.delete_sphere('SPHERE', self.form.cleft_partition_radius_slider))

        # NRGRank:
        self.form.nrgrank_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown_target(self.form.nrgrank_select_target, self.form.output_box))
        self.form.nrgrank_select_target.currentIndexChanged.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.nrgrank_select_binding_site, self.form.nrgrank_select_target.currentText(), self.form.output_box, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.nrgrank_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_select_ligand))
        self.form.nrgrank_button_start.clicked.connect(self.run_nrgrank)
        self.form.nrgrank_button_cancel.clicked.connect(self.abort_nrgrank)
        #TODO: in result table for FDA set show chemical name
        self.form.nrgrank_result_browse_button.clicked.connect(lambda: general_functions.folder_browser(self.form.nrgrank_result_path, os.path.join(self.form.temp_line_edit.text(), 'NRGRank'), "CSV file (*.csv)"))
        self.form.nrgrank_result_table.selectionModel().selectionChanged.connect(lambda: nrgrank_on_target.show_ligand_from_table(self.form.nrgrank_result_table, self.form.nrgrank_select_binding_site.currentText(), self.form.nrgrank_select_ligand.currentText()))

        # NRGRank ligand manager:
        self.form.nrgrank_add_ligandset_button.clicked.connect(lambda: general_functions.folder_browser(self.form.nrgrank_add_ligand_file_path, self.ligand_set_folder_path, 'All Files (*)'))
        self.form.nrgrank_delete_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_delete_ligand_set_dropdown, ignore_defaults=True))
        self.form.nrgrank_ligand_set_delete.clicked.connect(lambda: nrgrank_smiles_management.delete_ligand_set(self.form.nrgrank_delete_ligand_set_dropdown.currentText(), self.ligand_set_folder_path, self.form.output_box))
        self.form.nrgrank_button_ligandset_add.clicked.connect(self.run_generate_conformers)

        # FlexAID:
        self.form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown_target(self.form.flexaid_select_target,self.form.output_box))
        self.form.flexaid_select_target.currentIndexChanged.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.flexaid_select_binding_site, self.form.flexaid_select_target.currentText(), self.form.output_box, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.flexaid_select_ligand, self.form.output_box))
        self.form.flexaid_button_start.clicked.connect(self.run_flexaid)
        self.form.flexaid_button_pause.clicked.connect(lambda: pause_resume_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path))
        self.form.flexaid_button_stop.clicked.connect(lambda: stop_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path))
        self.form.flexaid_button_abort.clicked.connect(lambda: abort_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path, self.flexaid_manager))
        self.form.flexaid_result_table.selectionModel().selectionChanged.connect(lambda: flexaid_show_ligand_from_table(self.form))

        # Surfaces
        self.form.surfaces_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_object_1, self.form.output_box))
        self.form.surfaces_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_ligand_object_1, self.form.output_box, lig=1, add_none=1, prefer_none=True, no_warning=True))
        self.form.surfaces_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_object_2, self.form.output_box, add_none=1))
        self.form.surfaces_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_ligand_object_2, self.form.output_box, lig=1, add_none=1, prefer_none=True, no_warning=True))
        self.form.surfaces_button_run.clicked.connect(lambda: run_Surfaces.load_surfaces(self.form, self.form.temp_line_edit.text(), install_dir, self.binary_folder_path, self.binary_suffix))
        #issue with choosing which table to show see spike case study pt 2
        #TODO: Center ligand when ligand is specified and object 2 is empty
        self.form.surface_select_individual_result.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_individual_result.currentText() + '.txt')))
        self.form.surface_select_cf_comparison.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_cf_comparison.currentText() + '.csv')))
        self.form.surfaces_refresh_result.clicked.connect(lambda: run_Surfaces.refresh_res(self.form, os.path.join(self.form.temp_line_edit.text(), 'Surfaces')))
        self.form.surfaces_refresh_result.clicked.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_cf_comparison.currentText() + '.csv')))
        self.form.surfaces_button_interface.clicked.connect(lambda: run_Surfaces.read_and_select_residues(os.path.join(self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_individual_result.currentText() + '.txt'),
            self.form.surface_select_individual_result.currentText().split('_')[1], num_rows=self.form.TOPN_lineEdit_2.text()))

        # NRGTEN
        self.form.NRGten_target_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target_object_1, self.form.output_box))
        self.form.NRGten_target_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.NRGten_select_ligand_object_1, self.form.output_box, no_warning=True, lig=1, add_none=1, prefer_none=True))
        self.form.NRGten_target_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target_object_2, self.form.output_box, add_none=1, prefer_mutants=True))
        self.form.NRGten_dynasig_run.clicked.connect(self.run_NRGTen)
        self.form.NRGten_conf_ensem_run.clicked.connect(lambda: run_NRGTEN.conformational_ensemble(self.form, install_dir))

        # Single Mutations
        self.form.Modeller_refresh_object.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.Modeller_select_object, self.form.output_box))
        self.form.Modeller_refresh_residue.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.Modeller_select_residue, self.form.output_box, lig=1))
        self.form.modeller_checkbox_all.clicked.connect(lambda: run_modeller.check_all(self.form))
        self.form.modeller_button_mutate.clicked.connect(self.run_single_mutation)

        # IsoMIF
        self.form.ISOMIF_target_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target_object_1, self.form.output_box))
        self.form.ISOMIF_target_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target_object_2, self.form.output_box, add_none=1))
        self.form.ISOMIF_cleft_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.ISOMIF_select_cleft_object_1, self.form.ISOMIF_select_target_object_1.currentText(), self.form.output_box, show_all_objects=True))
        self.form.ISOMIF_cleft_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_ligand_object_1, self.form.output_box, lig=1, add_none=1))
        self.form.ISOMIF_cleft_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.ISOMIF_select_cleft_object_2, self.form.ISOMIF_select_target_object_2.currentText(), self.form.output_box, add_none=True, show_all_objects=True))
        self.form.ISOMIF_cleft_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_ligand_object_2, self.form.output_box, lig=1, add_none=1))
        self.form.ISOMIF_run.clicked.connect(lambda: run_isomif.mif_plot(self.form, self.binary_folder_path, self.binary_suffix, install_dir))

        # Settings functions
        self.form.Settings_button_refresh_obj.clicked.connect(lambda: run_settings.refresh_objects(self.form))
        self.form.Settings_button_combine_obj.clicked.connect(lambda: run_settings.combine_objects(self.form))
        self.form.Settings_split_states_button.clicked.connect(lambda: run_settings.devide_states(self.form))
        #TODO: Add option to remove all except selected
        self.form.Settings_refresh_hetatm.clicked.connect(lambda: run_settings.refresh_hetatm(self.form))
        self.form.Settings_button_remove_hetatm.clicked.connect(lambda: run_settings.remove_het_atoms(self.form))
        self.form.adv_ga_edit_button.clicked.connect(self.settings_edit_all_dialog)

    def settings_edit_all_dialog(self):
        from ga_settings_editor import SettingsGUI
        if self.advanced_settings_dialog is None:
            flexaid_ga_path = os.path.join(install_dir, 'deps', 'flexaid', 'ga_inp.dat')
            self.advanced_settings_dialog = SettingsGUI(flexaid_ga_path, self.form)
        self.advanced_settings_dialog.show()

    def run_getcleft(self):
        try:
            self.getcleftrunner = getcleft.GetCleftRunner(self.form, self.binary_folder_path, self.binary_suffix, install_dir, self.color_list)
        except ValueError as e:
            general_functions.output_message(self.form.output_box, e, 'warning')
        else:
            self.getcleftrunner.run_task()

    def run_nrgrank(self):
        self.nrgrankrunner = nrgrank_on_target.NRGRankManager(self.form, install_dir, self.ligand_set_folder_path, self.model)
        self.nrgrankrunner.run_nrgrank()

    def abort_nrgrank(self):
        self.nrgrankrunner.handle_thread_finished()
        self.nrgrankrunner = None

    def run_generate_conformers(self):
        smiles_path = self.form.nrgrank_add_ligand_file_path.text()
        if smiles_path == '':
            general_functions.output_message(self.form.output_box, 'Missing Smiles file path', 'warning')
            return
        self.conformer_generator = nrgrank_smiles_management.ConfGeneratorManager(self.form, install_dir, self.ligand_set_folder_path)
        self.conformer_generator.generate_conformer()

    def run_flexaid(self):
        self.flexaid_manager = FlexAIDManager(self.form, self.binary_folder_path, self.binary_suffix, install_dir, self.color_list, self.model)
        self.flexaid_manager.start_run()

    def run_NRGTen(self):
        general_functions.output_message(self.form.output_box, "=========== DynaSig ===========", 'valid')
        target_1 = self.form.NRGten_select_target_object_1.currentText()
        if target_1 == '':
            general_functions.output_message(self.form.output_box, 'No Object selected under: Load Object', 'warning')
        else:
            self.nrgtenrunner = DynasigManager(self.form, install_dir)
            self.nrgtenrunner.run_nrgten()

    def run_single_mutation(self):
        general_functions.output_message(self.form.output_box, '=========== Single Mutations ===========', 'valid')
        object_names = ['Object to mutate', 'Selected residue(s)']
        objects_to_check = ['Modeller_select_object', 'Modeller_select_residue']
        ok_continue = True
        for obj_counter, obj in enumerate(objects_to_check):
            obj_attribute = getattr(self.form, obj)
            if obj_attribute.currentText() == '':
                general_functions.output_message(self.form.output_box, f'{object_names[obj_counter]} not selected. Cannot run', 'warning')
                ok_continue = False
        if ok_continue:
            self.single_mutation_runner = single_mutationManager(self.form, install_dir)
            self.single_mutation_runner.run_single_mutation()


class NRGSuitePlugin(QWidget):
    def __init__(self):
        super().__init__()
        self.form = loadUi(os.path.join(install_dir, 'plugin.ui'), self)
        self.binary_suffix = None
        self.operating_system = None
        self.get_os()
        self.binary_folder_path = os.path.join(install_dir, 'bin', self.operating_system)
        test_binary(self.binary_folder_path, self.operating_system)
        self.get_folders()
        self.manage_dirs()
        self.check_modeller()
        self.set_flexaid_ga_params()
        self.form.stackedWidget.setCurrentIndex(0)
        self.form.flexaid_tab.setTabEnabled(2, False)
        self.form.NRGRank_tabs.setTabEnabled(2, False)
        general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box, no_warning=True)
        general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_select_ligand)
        self.form.nrgrank_cpu_usage_target.setCurrentText("75%")
        self.color_list = general_functions.load_color_list(os.path.join(install_dir, 'deps', 'getcleft', 'color_list.txt'))
        self.form.nrgrank_progress_label.setText('')
        self.form.nrgrank_loading_gif.setText('')
        self.form.nrgrank_progress.hide()
        self.form.mutation_progress_label.setText('')
        self.form.mutation_progress.hide()
        self.controller = Controller(self.form, self.binary_folder_path, self.binary_suffix, self.operating_system, self.ligand_set_folder_path, self.color_list)

    def get_os(self):
        operating_system = platform.system().upper()
        self.binary_suffix = ''
        if operating_system == 'LINUX' or operating_system == 'BSD':
            self.operating_system = 'linux'
        elif operating_system == 'DARWIN':
            self.operating_system = 'mac'
        elif operating_system == 'WINDOWS' or operating_system == 'MICROSOFT' or operating_system == 'WIN32':
            self.operating_system = 'win'
        else:
            exit('Unknown operating system')

    def get_folders(self):
        self.ligand_set_folder_path = os.path.join(install_dir, 'nrgrank_ligand_sets')
        self.plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
        self.temp_path = os.path.join(self.plugin_tmp_output_path, 'temp')
        self.form.temp_line_edit.setText(self.temp_path)
        self.nrgrank_output_path = os.path.join(self.form.temp_line_edit.text(), 'NRGRank')
        self.surfaces_output_path = os.path.join(self.form.temp_line_edit.text(), 'Surfaces')
        self.modeller_save_path = os.path.join(self.form.temp_line_edit.text(), 'modeller')
        self.nrgten_save_path = os.path.join(self.form.temp_line_edit.text(), 'NRGTEN')
        self.isomif_save_path = os.path.join(self.form.temp_line_edit.text(), 'IsoMIF')
        self.obj_manager_save_path = os.path.join(self.form.temp_line_edit.text(), 'Object_manager')

    def manage_dirs(self):
        if os.path.isdir(self.plugin_tmp_output_path):
            shutil.rmtree(self.plugin_tmp_output_path)
        os.mkdir(self.plugin_tmp_output_path)
        os.mkdir(self.form.temp_line_edit.text())
        os.mkdir(self.surfaces_output_path)
        os.mkdir(self.nrgrank_output_path)
        os.mkdir(self.modeller_save_path)
        os.mkdir(self.nrgten_save_path)
        os.mkdir(self.isomif_save_path)
        os.mkdir(self.obj_manager_save_path)

    def check_modeller(self):
        if 'modeller' not in sys.modules:
            general_functions.output_message(self.form.output_box, 'Modeller install not detected. '
                                                                   'The modeller tab will be unavailable. '
                                                                   'It will not possible to optimise states in the NRGTEN tab. '
                                                                   '\nPlease install via conda.', 'warning')
            general_functions.output_message(self.form.output_box, '=====================', 'warning')
            self.form.NRGten_optmizestates.setEnabled(False)
            self.form.button_modeller.setEnabled(False)
            self.form.button_modeller.setStyleSheet("background-color: black; color: white;")

    def set_flexaid_ga_params(self):
        flexaid_ga_path = os.path.join(install_dir, 'deps', 'flexaid', 'ga_inp.dat')
        with open(flexaid_ga_path, 'r') as f:
            for line in f:
                if line.startswith('NUMCHROM'):
                    value = line.strip().split()[-1]
                    self.form.input_num_chromosomes.setText(str(value))
                if line.startswith('NUMGENER'):
                    value = line.strip().split()[-1]
                    self.form.input_num_generations.setText(str(value))