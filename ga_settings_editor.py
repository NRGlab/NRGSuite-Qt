import os.path
import sys, shutil
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QTableWidget,
    QTableWidgetItem, QLabel
)

class SettingsGUI(QWidget):
    def __init__(self, table_path, form):
        super().__init__()
        self.setWindowTitle("GA Settings Editor")
        self.resize(300, 400)
        self.table_path = table_path
        self.form = form
        self.default_table_path = os.path.join(os.path.dirname(self.table_path), 'ga_inp_default.dat')
        self.layout = QVBoxLayout()
        self.table = QTableWidget(0, 2)
        self.layout.addWidget(self.table)
        self.load_settings()
        self.formatTable()

        self.default_button = QPushButton("Reset Default")
        self.default_button.clicked.connect(self.reset_to_default)
        self.layout.addWidget(self.default_button)

        self.layout.addWidget(QLabel("Settings are saved when closing this window"))
        self.setLayout(self.layout)
        self.table.itemChanged.connect(self.save_settings)

    def load_settings(self):
        with open(self.table_path, "r") as file:
            lines = file.readlines()
        self.table.setRowCount(len(lines))
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if parts:
                self.table.setItem(i, 0, QTableWidgetItem(parts[0]))
                self.table.setItem(i, 1, QTableWidgetItem(" ".join(parts[1:])))

    def formatTable(self):
        hh = self.table.horizontalHeader()
        hh.setStretchLastSection(True)
        hh.setVisible(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setFocus()
        self.table.hide()
        self.table.resizeColumnsToContents()
        self.table.show()

    def save_settings(self):
        with open(self.table_path, "w") as file:
            for row in range(self.table.rowCount()):
                param = self.table.item(row, 0).text() if self.table.item(row, 0) else ""
                value = self.table.item(row, 1).text() if self.table.item(row, 1) else ""
                file.write(f"{param} {value}\n")
                if param == 'NUMCHROM':
                    self.form.input_num_chromosomes.setText(str(value))
                if param == 'NUMGENER':
                    self.form.input_num_generations.setText(str(value))

    def reset_to_default(self):
        os.remove(self.table_path)
        shutil.copy(self.default_table_path, self.table_path)
        self.load_settings()
        self.formatTable()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SettingsGUI('/opt/miniconda3/lib/python3.12/site-packages/pmg_tk/startup/NRGSuite-Qt/deps/flexaid/ga_inp.dat', '')
    window.show()
    sys.exit(app.exec())
