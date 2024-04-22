import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout
from PyQt5.QtCore import QPropertyAnimation, QUrl
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.uic import loadUi
from conect import get_info
from PyQt5.QtWebEngineWidgets import QWebEngineView
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import os
import threading
import http.server


class CORSHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()

def start_server():
    httpd = http.server.HTTPServer(("", 8001), CORSHTTPRequestHandler)
    httpd.serve_forever()

threading.Thread(target=start_server, daemon=True).start()


class VentanaPrincipal(QMainWindow):
    def __init__(self):
        super(VentanaPrincipal, self).__init__()
        loadUi('phylo.ui', self)

        self.btn_menu.clicked.connect(self.slide_menu)
        self.btn_restaurar.clicked.connect(self.restaurar)
        self.btn_close.clicked.connect(lambda: self.close())

        self.btn_menu.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.btn_restaurar.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.btn_close.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        self.NavBar.mouseMoveEvent = self.mover_ventana

        #elimina barra de titulo -opacidad
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        self.setWindowOpacity(1)

        self.gripSize = 10
        self.grip = QtWidgets.QSizeGrip(self)
        self.grip.resize(self.gripSize, self.gripSize)

        self.datos = get_info()
        # print(datos)

        self.ListGen.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

        # Crea el modelo para la QListView.
        model = QtGui.QStandardItemModel(self.ListGen)

        # Añade los datos al modelo.
        for item in self.datos:
            model.appendRow(QtGui.QStandardItem(item['organism']))

        
        # Asigna el modelo a la QListView.
        self.ListGen.setModel(model)

        self.ListGen.clicked.connect(self.on_item_clicked)

        
        '''Canvas Widget Config'''
        self.view = QWebEngineView()
        layout = QVBoxLayout()
        self.Visor.setLayout(layout)
        layout.addWidget(self.view)

    def on_item_clicked(self, index):
        # Obtiene los datos del elemento seleccionado.
        item_data = self.datos[index.row()]

        # Llama a plot con los datos del elemento seleccionado.
        self.plot(item_data)

    def plot(self, data):
        # Procesa los datos.
        dna_sequence = data['info']

        # Convierte la secuencia de ADN en una secuencia de aminoácidos
        protein_seq = Seq(dna_sequence).translate(to_stop=False)
        protein_seq = protein_seq.strip('*')

        # Convierte la secuencia de aminoácidos en una secuencia de tres letras
        three_letter_seq = seq3(protein_seq)

        # Crea un nuevo objeto Structure
        structure = Structure.Structure("My_Protein")

        # Añade un modelo al objeto Structure
        model = Model.Model(0)
        structure.add(model)

        # Añade una cadena al modelo
        chain = Chain.Chain("A")
        model.add(chain)

        # Añade los residuos a la cadena
        for i, residue in enumerate(three_letter_seq.split("-")):
            res = Residue.Residue((" ", i, " "), residue, " ")
            chain.add(res)

        # Guarda el objeto Structure en un archivo
        io = PDBIO()
        io.set_structure(structure)
        pdb_file_name = "my_protein.pdb"
        io.save(pdb_file_name)

        # Actualiza la vista del visor con el nuevo archivo PDB
        self.update_viewer(pdb_file_name)

    def update_viewer(self, pdb_file_name):
        # Carga el archivo HTML en tu aplicación PyQt5
        self.view.load(QUrl(f"http://localhost:8001/index.html?pdb={pdb_file_name}"))




    def restaurar(self):
        if self.isMaximized() == False:
            self.showMaximized()
            self.btn_restaurar.setText('Min')
        else:
            self.showNormal()
            self.btn_restaurar.setText('Max')

    #mover ventana
    def mousePressEvent(self, event):
        self.click_position = event.globalPos()

    def mover_ventana(self, event):
        if self.isMaximized() == False:
            if event.buttons()== QtCore.Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.click_position)
                self.click_position = event.globalPos()
                event.accept()

        if event.globalPos().y() <=10:
            self.showMaximized()
        else:
            self.showNormal()
        
    def slide_menu(self):
        if True:
            width = self.Controls.width()
            normal = 0
            if width == 0:
                extender = 230
            else:
                extender = normal
            self.animacion = QPropertyAnimation(self.Controls, b'minimumWidth')
            self.animacion.setDuration(300)
            self.animacion.setStartValue(width)
            self.animacion.setEndValue(extender)
            self.animacion.setEasingCurve(QtCore.QEasingCurve.InOutQuart)
            self.animacion.start()




if __name__ == "__main__":
    app = QApplication(sys.argv)
    mi_app = VentanaPrincipal()
    mi_app.show()
    sys.exit(app.exec_())