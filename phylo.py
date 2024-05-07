from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout
from PyQt5.QtCore import QPropertyAnimation
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import numpy as np
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

class CustomToolbar(NavigationToolbar2QT):
    toolitems = [t for t in NavigationToolbar2QT.toolitems if t[0] == 'Save']

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
        self.ax = self.fig.add_subplot(111, projection='3d')

        # Agregamos el gráfico a la ventana
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.canvas.updateGeometry()

        layout.addWidget(self.canvas)
        self.toolbar = CustomToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)

        # Tu secuencia de ADN
        self.my_seq = Seq(dna_sequence)

        # Dibujamos las líneas entre las bases y sus complementos
        self.draw_dna()

        # Conectamos el evento scroll_event a la función on_scroll
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    def on_scroll(self, event):
        # Zoom in or out
        base_scale = 0.5
        # get the current x and y limits
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        xdata = event.xdata  # get event x location
        ydata = event.ydata  # get event y location
        if event.button == 'up':
            # deal with zoom out
            scale_factor = base_scale
        elif event.button == 'down':
            # deal with zoom in
            scale_factor = 1/base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
        # set new limits
        self.ax.set_xlim([xdata - cur_xrange*scale_factor,
                        xdata + cur_xrange*scale_factor])
        self.ax.set_ylim([ydata - cur_yrange*scale_factor,
                        ydata + cur_yrange*scale_factor])
        self.canvas.draw()



    def draw_dna(self):
        # Limpiamos el gráfico
        self.ax.clear()
        self.ax.axis('off')

        # Establecemos el color de fondo del gráfico a negro
        self.ax.set_facecolor('black')

        # Obtenemos la secuencia complementaria
        complement_seq = self.my_seq.complement()

        # Verificamos si la secuencia es válida
        if not all(base in ["A", "T", "G", "C"] for base in self.my_seq):
            text = 'This sequence is not valid'
            self.ax.text(0.5, 0.5, 0.5, text, va='center', ha='center', color='#fff')
            return

        # Dibujamos las líneas curvas entre las bases y sus complementos
        seq_len = len(self.my_seq)
        for pos, (base, comp_base) in enumerate(zip(self.my_seq, complement_seq)):
            x = np.linspace(pos+1, pos+1, 100)
            y_base = 3 + 0.02 * min(pos, seq_len - pos)
            y_comp = 6 - 0.02 * min(pos, seq_len - pos)
            y = np.linspace(y_base, y_comp, 100)
            z = 0.1 * (x-pos-1)**2 + pos+1
            color = ''
            if base == "A":
                color = '#ff3535ff'
            elif base == "T":
                color = '#00ff57ff'
            elif base == "G":
                color = '#ffe600ff'
            elif base == "C":
                color = '#00efffff'
            self.ax.plot(x[:50], y[:50], z[:50], color=color)
            if comp_base == "A":
                color = '#ff3535ff'
            elif comp_base == "T":
                color = '#00ff57ff'
            elif comp_base == "G":
                color = '#ffe600ff'
            elif comp_base == "C":
                color = '#00efffff'
            self.ax.plot(x[50:], y[:50], z[:50:], color=color)

        # Dibujamos las líneas que conectan las bases de la misma secuencia
        x = np.arange(1, seq_len + 1)
        y = 3 + 0.02 * np.minimum(x, seq_len - x + 1)
        z = 0.1 * (x-x)**2 + x
        self.ax.plot(x, y, z, color='white')
        y = 6 - 0.02 * np.minimum(x, seq_len - x + 1)
        self.ax.plot(x, y, z, color='white')

        # Ajustamos los límites del gráfico
        self.ax.set_zlim(0, len(self.my_seq)+1)
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
        self.ax.axis('off')

        # Establecemos el color de fondo del gráfico a negro
        self.ax.set_facecolor('black')

        # Obtenemos la secuencia complementaria
        complement_seq = self.my_seq.complement()

        # Verificamos si la secuencia es válida
        if not all(base in ["A", "T", "G", "C"] for base in self.my_seq):
            text = 'This sequence is not valid'
            self.ax.text(0.5, 0.5, 0.5, text, va='center', ha='center', color='#fff')
            return

        # Dibujamos las líneas curvas entre las bases y sus complementos
        seq_len = len(self.my_seq)
        for pos, (base, comp_base) in enumerate(zip(self.my_seq, complement_seq)):
            x = np.linspace(pos+1, pos+1, 100)
            y_base = 3 + 0.02 * min(pos, seq_len - pos)  # Las bases están en un orden escalonado en el eje y
            y_comp = 6 - 0.02 * min(pos, seq_len - pos)  # Los complementos están en un orden escalonado en el eje y
            y = np.linspace(y_base, y_comp, 100)
            z = 0.1 * (x-pos-1)**2 + pos+1
            color = ''
            if base == "A":
                color = '#ff3535ff'
            elif base == "T":
                color = '#00ff57ff'
            elif base == "G":
                color = '#ffe600ff'
            elif base == "C":
                color = '#00efffff'
            self.ax.plot(x[:50], y[:50], z[:50], color=color)
            if comp_base == "A":
                color = '#ff3535ff'
            elif comp_base == "T":
                color = '#00ff57ff'
            elif comp_base == "G":
                color = '#ffe600ff'
            elif comp_base == "C":
                color = '#00efffff'
            self.ax.plot(x[50:], y[50:], z[50:], color=color)

        # Dibujamos las líneas que conectan las bases de la misma secuencia
        x = np.arange(1, seq_len + 1)
        y = 3 + 0.02 * np.minimum(x, seq_len - x + 1)
        z = 0.1 * (x-x)**2 + x
        self.ax.plot(x, y, z, color='white')
        y = 6 - 0.02 * np.minimum(x, seq_len - x + 1)
        self.ax.plot(x, y, z, color='white')

        # Ajustamos los límites del gráfico
        self.ax.set_zlim(0, len(self.my_seq)+1)
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