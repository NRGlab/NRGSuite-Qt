from PyQt5.QtCore import pyqtSignal, QThread, Qt
import general_functions
import subprocess
import os
from PyQt5.QtGui import QColor
from pymol import cmd
import datetime
import time
from PyQt5.QtGui import QStandardItem

def pause_simulation(form, specific_simulation_path):
    with open(os.path.join(specific_simulation_path, '.pause'), 'a'):
        pass
    form.flexaid_button_pause.setText('Resume')

def resume_simulation(form, specific_simulation_path):
    os.remove(os.path.join(specific_simulation_path, '.pause'))
    form.flexaid_button_pause.setText('Pause')

def abort_simulation(form, specific_simulation_path, flexaid_manager):
    with open(os.path.join(specific_simulation_path, '.abort'), 'a'):
        pass
    form.flexaid_button_start.setEnabled(True)
    form.flexaid_tab.setCurrentIndex(1)
    flexaid_manager.quit_workers()

def stop_simulation(form, specific_simulation_path):
    with open(os.path.join(specific_simulation_path, '.stop'), 'a'):
        pass
    if form.flexaid_button_pause.text() == 'Resume':
        resume_simulation(form, specific_simulation_path)
    form.flexaid_button_start.setEnabled(True)

def pause_resume_simulation(form, specific_simulation_path):
    if form.flexaid_button_pause.text() == 'Pause':
        pause_simulation(form, specific_simulation_path)
    elif form.flexaid_button_pause.text() == 'Resume':
        resume_simulation(form, specific_simulation_path)


def flexaid_show_ligand_from_table(form):
    table_object = form.flexaid_result_table
    binding_site = form.flexaid_select_binding_site.currentText()
    ligand = form.flexaid_select_ligand.currentText()

    selected_indexes = table_object.selectionModel().selectedIndexes()
    ligands_to_show = []
    setting_dictionary = {'number_chromosomes': int(form.input_num_chromosomes.text()),
                          'number_generations': int(form.input_num_generations.text()),
                          'max_results': int(form.flexaid_max_results.text())}
    for index in selected_indexes:
        cell_text = index.data()
        column = index.column()
        if column == 1:
            ligands_to_show.append(str(int(cell_text) - 1))
    if ligands_to_show:
        if ligand.find('bd_site') != -1:
            ligand = ligand.split('_')[0]
        ligands_to_show = [f"RESULT_{item}_flx_{ligand}_{binding_site}_{setting_dictionary['number_chromosomes']}x{setting_dictionary['number_generations']}" for item in ligands_to_show]
        group_name = general_functions.get_group_of_object(ligands_to_show[0])
        if not group_name:
            print('no group found')
            return
        else:
            group_objects = cmd.get_object_list(f'({group_name})')
            for group_object in group_objects:
                if group_object not in ligands_to_show:
                    cmd.disable(group_object)
                else:
                    cmd.enable(group_object)


class FlexAIDManager:
    def __init__(self, form, binary_folder_path, binary_suffix, install_dir, color_list, model):
        self.form = form
        self.model = model
        self.flexaid_temp_path = os.path.join(self.form.temp_line_edit.text(), 'FlexAID')
        self.flexaid_simulation_folder = os.path.join(self.flexaid_temp_path, 'Simulation')
        os.makedirs(self.flexaid_temp_path, exist_ok=True)
        os.makedirs(self.flexaid_simulation_folder, exist_ok=True)
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.install_dir = install_dir
        self.rmsd = form.flexaid_ligref_checkBox.isChecked()
        self.run_specific_simulate_folder_path = os.path.join(self.flexaid_simulation_folder, datetime.datetime.now().strftime("%d-%m-%y-%I-%M-%S"))
        self.update_file_path = os.path.join(self.run_specific_simulate_folder_path, ".update")
        self.simulation_settings = self.get_simulation_settings()
        self.max_generations = int(self.simulation_settings['number_generations'])
        self.form.flexaid_progress.setValue(0)
        self.form.flexaid_progress.setMaximum(self.max_generations)
        self.form.generation_label.setText(f'Generation: 0/{self.max_generations}')
        self.flexaid_deps_path = os.path.join(install_dir, 'deps', 'flexaid')
        self.color_list = color_list
        self.model.clear()
        self.target_name = self.form.flexaid_select_target.currentText()
        self.target_save_path = os.path.join(self.flexaid_temp_path, 'flexaid_target.pdb')
        self.threads = []

    def start_run(self):
        self.start_file_updater()
        self.start_flexaid_thread()

    def start_file_updater(self):
        self.file_updater_thread = FileUpdaterThread(self.update_file_path)
        self.file_updater_thread.update_signal.connect(self.perform_updates)
        self.threads.append(self.file_updater_thread)
        self.file_updater_thread.start()

    def perform_updates(self, current_generation, data):
        self.form.flexaid_progress.setValue(current_generation)
        self.form.generation_label.setText(f'Generation: {current_generation}/{self.max_generations}')
        self.colour_specific_cell(data)

    def start_flexaid_thread(self):
        cmd.save(self.target_save_path, self.target_name)
        self.ligand_name = self.form.flexaid_select_ligand.currentText()
        ligand_save_path = os.path.join(self.flexaid_temp_path, 'flexaid_ligand.pdb')
        cmd.save(ligand_save_path, self.ligand_name)
        self.binding_site_name = self.form.flexaid_select_binding_site.currentText()
        binding_site_path = os.path.join(self.flexaid_temp_path, 'binding_site_sph_.pdb')
        cmd.save(binding_site_path, self.binding_site_name)
        self.toggle_buttons(True)
        general_functions.output_message(self.form.output_box, f"=========== FlexAID ===========", 'valid')
        general_functions.output_message(self.form.output_box, 'Docking...', 'valid')
        self.form.flexaid_tab.setCurrentIndex(2)
        self.form.flexaid_tab.setTabEnabled(2, True)
        self.worker = FlexAIDWorkerThread(self.binary_folder_path, self.binary_suffix, self.install_dir,
                                          self.run_specific_simulate_folder_path, self.flexaid_temp_path,
                                          self.target_save_path, ligand_save_path, binding_site_path, self.simulation_settings)
        self.worker.finished.connect(self.flexaid_end_processes)
        self.threads.append(self.worker)
        self.worker.start()

    def get_simulation_settings(self):
        setting_dictionary = {'number_chromosomes': int(self.form.input_num_chromosomes.text()),
                              'number_generations': int(self.form.input_num_generations.text()),
                              'max_results': int(self.form.flexaid_max_results.text())}
        return setting_dictionary


    def colour_specific_cell(self, data):
        table_headers = ['Colour', 'TOP', 'CF', 'Fitness', 'Last RMSD']
        self.model.clear()
        self.model.setHorizontalHeaderLabels(table_headers)
        for row in data:
            item_list = []
            for header in table_headers:
                if header == 'Colour':
                    item = QStandardItem()
                    r, g, b = self.color_list[row['TOP'] - 1]["rgb"]
                    item.setBackground(QColor(r, g, b))
                else:
                    row_item = row[header]
                    if type(row_item) == float:
                        row_item = f"{row_item:.2f}"
                    else:
                        row_item = str(row_item)
                    item = QStandardItem(row_item)
                item.setTextAlignment(Qt.AlignCenter)
                item_list.append(item)
            self.model.appendRow(item_list)
        self.form.flexaid_result_table.setModel(self.model)
        self.form.flexaid_result_table.resizeColumnsToContents()
        padding = 20
        for column in range(self.form.flexaid_result_table.model().columnCount()):
            current_width = self.form.flexaid_result_table.columnWidth(column)
            self.form.flexaid_result_table.setColumnWidth(column, current_width + padding)


    def final_table_update(self):
        try:
            with open(os.path.join(os.path.join(self.run_specific_simulate_folder_path, 'result.cad'))) as f:
                lines = f.readlines()
        except FileNotFoundError:
            general_functions.output_message(self.form.output_box, f"No result file: {os.path.join(self.run_specific_simulate_folder_path, 'result.cad')}", 'error')
        result_counter = 1
        data = []
        for line in lines:
            if line.startswith('Cluster'):
                line = line.split(':')[-1].strip().split(' ')
                cf = float(line[1].split('=')[-1])
                rmsd = 'N/A'
                if self.rmsd:
                    with open(os.path.join(self.run_specific_simulate_folder_path, f'RESULT_{result_counter-1}.pdb'), 'r') as result_file:
                        lines = result_file.readlines()
                        for line in lines:
                            if 'RMSD to' in line:
                                rmsd = float(line.split()[1])
                data.append({'Colour': result_counter - 1, 'TOP': result_counter, 'CF': cf, 'Fitness': 'N/A', 'Last RMSD': rmsd})
                result_counter += 1
        self.colour_specific_cell(data)


    def quit_workers(self):
        for worker in self.threads:
            worker.stop()
            worker.quit()
            worker.wait()


    def flexaid_end_processes(self):
        general_functions.output_message(self.form.output_box, f"=========== FlexAID ===========", 'valid')
        self.quit_workers()
        self.toggle_buttons(False)
        self.final_table_update()
        self.load_show_flexaid_result()
        self.form.flexaid_progress.setValue(self.max_generations)
        self.form.generation_label.setText(f'Generation: {self.max_generations}/{self.max_generations}')

    def toggle_buttons(self, true_false):
        self.form.flexaid_button_pause.setEnabled(true_false)
        self.form.flexaid_button_abort.setEnabled(true_false)
        self.form.flexaid_button_stop.setEnabled(true_false)

    def load_show_flexaid_result(self):
        result_files = next(os.walk(self.run_specific_simulate_folder_path), (None, None, []))[2]
        result_files = [file for file in result_files if file.startswith('RESULT') and not file.endswith('INI.pdb') and file.endswith('.pdb')]
        result_files = sorted(result_files, key=lambda s: int(s.split('_')[1].split('.')[0]))
        cmd.disable('everything')
        true_result = 0
        cmd.delete(f'FlexAID_{self.binding_site_name}')
        cmd.disable(f'{self.target_name},NRGRank,GetCleft,Surfaces')
        if self.ligand_name.find('bd_site') != -1:
            self.ligand_name = self.ligand_name.split('_')[0]
        results_group_name = f"flx_{self.ligand_name}_{self.binding_site_name}_{self.simulation_settings['number_chromosomes']}x{self.simulation_settings['number_generations']}"
        for file in result_files:
            file_path = os.path.join(self.run_specific_simulate_folder_path, file)
            molecule_name = f"{os.path.basename(file_path[:-4])}_{results_group_name}"
            cmd.load(file_path, molecule_name)
            color_name = self.color_list[true_result]['name']
            cmd.color('white', molecule_name)
            cmd.color(color_name, f"{molecule_name} and resn LIG")
            cmd.group(results_group_name, molecule_name)
            true_result += 1
        cmd.group('FlexAID', results_group_name)


class FileUpdaterThread(QThread):
    update_signal = pyqtSignal(int, list)

    def __init__(self, update_file_path):
        super(FileUpdaterThread, self).__init__()
        self.update_file_path = update_file_path
        self.current_generation = 0
        self._is_running = True

    def stop(self):
        self._is_running = False

    def run(self):
        while self._is_running:
            if os.path.exists(self.update_file_path):
                self.read_update()
                self.current_generation += 1
            time.sleep(0.01)

    def read_update(self):
        data = []
        if self.current_generation % 10 == 0:
            with open(self.update_file_path, "r") as f:
                for line_counter, line in enumerate(f):
                    if line_counter > 1:
                        line = line.split()
                        top_number = int(line[0]) + 1
                        cf = ''
                        for item in line:
                            if item.startswith('cf='):
                                cf = item.split('=')[-1]
                                if cf == '':
                                    cf = line[-5]
                                break
                        cf = float(cf)
                        fitness = float(line[-1])
                        rmsd = 'N/A'
                        data.append({'Colour': top_number - 1, 'TOP': top_number, 'CF': cf, 'Fitness': fitness, 'Last RMSD': rmsd})
            self.update_signal.emit(self.current_generation, data)
        for _ in range(5):
            try:
                os.remove(self.update_file_path)
                break
            except PermissionError:
                time.sleep(0.1)


class FlexAIDWorkerThread(QThread):
    finished = pyqtSignal()

    def __init__(self, binary_folder_path, binary_suffix, install_dir, flexaid_result_path, flexaid_temp_path,
                 target_save_path, ligand_save_path, binding_site_path, setting_dictionary):
        super(FlexAIDWorkerThread, self).__init__()
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.install_dir = install_dir
        self.flexaid_result_path = flexaid_result_path
        self.flexaid_temp_path = flexaid_temp_path
        self.target_save_path = target_save_path
        self.ligand_save_path = ligand_save_path
        self.binding_site_path = binding_site_path
        self.setting_dictionary = setting_dictionary
        self.max_results = setting_dictionary['max_results']
        self._is_running = True
        self._process = None  # Store the process handle

    @staticmethod
    def process_ligand(process_ligand_path, input_path, istarget=False):
        if istarget:
            process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" -target'
        else:
            process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" --atom_index 90000 -ref'
        subprocess.run(process_ligand_command, shell=True)

    @staticmethod
    def count_flex(ligand_inp_file_path):
        with open(ligand_inp_file_path, 'r') as t:
            texto = t.readlines()
            count = 0
            for line in texto:
                if line[0] == 'F':
                    count += 1
        return count

    def write_config(self, target_inp_path, cleft, ligand_inp_path, max_results, flexaid_output_path,
                     flexaid_result_path, flexaid_deps_path):
        with open(os.path.join(flexaid_deps_path, 'template_config.inp'), "r") as t1:
            lines = t1.readlines()
        config_file_output_path = os.path.join(flexaid_output_path, 'config.inp')
        matrix_path = os.path.join(flexaid_deps_path, 'MC_st0r5.2_6.dat')
        with open(config_file_output_path, "w") as t2:
            for line in lines:
                if line.startswith('PDBNAM'):
                    t2.write(f'PDBNAM {target_inp_path}\n')
                elif line.startswith('INPLIG'):
                    t2.write(f'INPLIG {ligand_inp_path}\n')
                elif line.startswith('RNGOPT LOCCLF'):
                    t2.write(f"RNGOPT LOCCLF {cleft}\n")
                elif line.startswith('OPTIMZ 9999 - 0'):
                    t2.write(line)
                    for flex in range(self.count_flex(ligand_inp_path)):
                        t2.write(f'OPTIMZ 9999 - {str(flex + 1)}\n')
                elif line.startswith('STATEP'):
                    t2.write(f'STATEP {flexaid_result_path}\n')
                elif line.startswith('RMSDST'):
                    t2.write(f'RMSDST {ligand_inp_path[:-4]}_ref.pdb\n')
                elif line.startswith('IMATRX'):
                    t2.write(f'IMATRX {matrix_path}\n')
                elif line.startswith('DEPSPA'):
                    t2.write(f'DEPSPA {flexaid_deps_path}\n')
                elif line.startswith('MAXRES'):
                    t2.write(f'MAXRES {max_results}\n')
                elif line.startswith('TEMPOP'):
                    t2.write(f'TEMPOP {os.path.join(flexaid_result_path, "temp")}\n')
                else:
                    t2.write(line)
        return config_file_output_path

    @staticmethod
    def edit_ga(ga_template_path, ga_write_path, setting_dictionary):
        with open(ga_template_path, "r") as ga_template:
            lines = ga_template.readlines()
        with open(ga_write_path, "w") as ga_write:
            for line in lines:
                if line.startswith('NUMCHROM'):
                    ga_write.write(f'NUMCHROM {setting_dictionary["number_chromosomes"]}\n')
                elif line.startswith('NUMGENER'):
                    ga_write.write(f'NUMGENER {setting_dictionary["number_generations"]}\n')
                else:
                    ga_write.write(line)

    def stop(self):
        self._is_running = False
        if self._process and self._process.poll() is None:
            self._process.terminate()

    def run(self):
        if self._is_running:
            flexaid_binary_path = os.path.join(self.binary_folder_path, f'FlexAID{self.binary_suffix}')
            process_ligand_path = os.path.join(self.binary_folder_path, f'Process_ligand{self.binary_suffix}')
            flexaid_deps_path = os.path.join(self.install_dir, 'deps', 'flexaid')
            os.mkdir(self.flexaid_result_path)
            os.mkdir(os.path.join(self.flexaid_result_path, 'temp'))
            flexaid_result_name_path = os.path.join(self.flexaid_result_path, "RESULT").replace('\\', '/')

            self.process_ligand(process_ligand_path, self.target_save_path, istarget=True)
            self.process_ligand(process_ligand_path, self.ligand_save_path)
            target_inp_path = os.path.splitext(self.target_save_path)[0] + '.inp.pdb'
            ligand_inp_path = os.path.splitext(self.ligand_save_path)[0] + '.inp'
            config_file_path = self.write_config(target_inp_path, self.binding_site_path, ligand_inp_path, self.max_results,
                                                 self.flexaid_temp_path, self.flexaid_result_path,
                                                 flexaid_deps_path).replace('\\', '/')

            ga_path = os.path.join(self.flexaid_temp_path, 'ga_inp.dat').replace('\\', '/')
            self.edit_ga(os.path.join(flexaid_deps_path, 'ga_inp.dat'), ga_path, self.setting_dictionary)
            flexaid_command = [flexaid_binary_path, config_file_path, ga_path, flexaid_result_name_path]

            self._process = subprocess.Popen(flexaid_command)
            while self._is_running and self._process.poll() is None:
                self.msleep(100)

            if self._process.poll() == 0:
                self.finished.emit()
