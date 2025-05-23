import os
import sys
import subprocess
import shutil
import importlib.metadata
import platform


def __init_plugin__(app):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('NRGSuite-Qt', run_plugin_gui)


def install_package(package, main_folder_path):
    try:
        __import__(package.split('=')[0])
    except ImportError as e:
        if package == 'modeller':
            print('MODELLER install not detected. Please install via conda. The MODELLER and NRGTen tabs will be unavailable')
        else:
            if package == 'Bio':
                package = 'biopython'
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable.replace('pythonw', 'python'), "-m", "pip", "install", package])
            if package == 'nrgten':
                distribution = importlib.metadata.distribution(package)
                shutil.copy(os.path.join(main_folder_path, 'deps', 'nrgten', 'amino_acids.atomtypes'),
                            os.path.join(str(distribution.locate_file('')), 'nrgten', 'config', 'amino_acids.atomtypes'))
                shutil.copy(os.path.join(main_folder_path, 'deps', 'nrgten', 'amino_acids.masses'),
                            os.path.join(str(distribution.locate_file('')), 'nrgten', 'config', 'amino_acids.masses'))


def check_packages(install_dir):
    print('Checking python modules')
    packages = ['nrgten', 'Bio', 'pandas', 'matplotlib', 'colour', 'scipy', 'numpy==2.0', 'numba', 'plotly', 'rdkit', 'joblib']
    for package in packages:
        install_package(package, install_dir)


def run_plugin_gui():
    install_dir = os.path.dirname(__file__)
    sys.path.append(install_dir)

    global dialog
    if dialog is None:
        check_packages(install_dir)
        import gui_main
        dialog = gui_main.NRGSuitePlugin()
        
    distribution = importlib.metadata.distribution('nrgten')
    vcon_exe_path = os.path.join(str(distribution.locate_file('')), 'pyvcon', 'vcontacts', 'vcon')
    if platform.system() == 'Windows' and (not os.path.isfile(vcon_exe_path) or os.path.getsize(vcon_exe_path)) != os.path.getsize(os.path.join(install_dir, 'bin', 'win', 'vcon_nrgten.exe')):
        shutil.copy(os.path.join(install_dir, 'bin', 'win', 'vcon_nrgten.exe'),
                    os.path.join(str(distribution.locate_file('')),'pyvcon', 'vcontacts', 'vcon'))
        if os.path.isfile(os.path.join(str(distribution.locate_file('')),'pyvcon', 'vcontacts', 'vcontacts.exe')):
            os.remove(os.path.join(str(distribution.locate_file('')),'pyvcon', 'vcontacts', 'vcontacts.exe'))


    dialog.show()
    dialog.raise_()

dialog = None

if __name__ == '__main__':
    from PyQt6.QtWidgets import QApplication
    import gui_main
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = gui_main.NRGSuitePlugin()
    window.show()
    sys.exit(app.exec_())
