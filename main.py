# En tu otro archivo
from phylo import Viewer
from PyQt5.QtWidgets import QApplication
import sys


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mi_app = Viewer("GCCAATGCAG")
    mi_app.show()
    sys.exit(app.exec_())