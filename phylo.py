from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout
from PyQt5.QtCore import QPropertyAnimation
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from Bio.Seq import Seq
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

class Viewer(QMainWindow):
    def __init__(self, dna_sequence):
        super(Viewer, self).__init__()
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
        
        layout = QVBoxLayout()
        self.Visor.setLayout(layout)

        # Creamos el gráfico
        self.fig = Figure(figsize=(10, 5), facecolor='black')
        self.ax = self.fig.add_subplot(111)

        # Agregamos el gráfico a la ventana
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setStyleSheet("background-color:black;")
        layout.addWidget(self.canvas)

        # Tu secuencia de ADN
        self.my_seq = Seq(dna_sequence)

        # Dibujamos las líneas entre las bases y sus complementos
        self.draw_dna()


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

    def draw_dna(self):
        # Limpiamos el gráfico
        self.ax.clear()

        # Obtenemos la secuencia complementaria
        complement_seq = self.my_seq.complement()

        # Dibujamos las líneas entre las bases y sus complementos
        for pos, (base, comp_base) in enumerate(zip(self.my_seq, complement_seq)):
            self.ax.plot([pos+1, pos+1], [3, 6], color='#ff00b8ff')
            if base == "A":
                self.ax.text(pos+1, 3, base, ha='center', va='center', color='#ff3535ff')
                self.ax.text(pos+1, 6, comp_base, ha='center', va='center', color='#00ff57ff')
            elif base == "T":
                self.ax.text(pos+1, 3, base, ha='center', va='center', color='#00ff57ff')
                self.ax.text(pos+1, 6, comp_base, ha='center', va='center', color='#ff3535ff')
            elif base == "G":
                self.ax.text(pos+1, 3, base, ha='center', va='center', color='#ffe600ff')
                self.ax.text(pos+1, 6, comp_base, ha='center', va='center', color='#00efffff')
            elif base == "C":
                self.ax.text(pos+1, 3, base, ha='center', va='center', color='#00efffff')
                self.ax.text(pos+1, 6, comp_base, ha='center', va='center', color='#ffe600ff')

        # Ajustamos los límites del gráfico
        self.ax.set_ylim(0, 9)
        self.ax.set_xlim(0, len(self.my_seq)+1)

        # Ajustamos las etiquetas del eje y para mostrar las bases en lugar de números
        self.ax.set_yticks([])
        self.ax.set_xticks([])

        # Cambiamos el color de fondo del gráfico a negro
        self.ax.set_facecolor('black')

        # Ocultamos los ejes
        self.ax.axis('off')

        # Redibujamos el gráfico
        self.canvas.draw()