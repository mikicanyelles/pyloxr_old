import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import uic
#import ctypes


class Menu_prin(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        uic.loadUi("ui/mainwindow.ui")        #Carrega el fitxer .ui

        self.setWindowTitle('Menu - pyLOXr')    #Canvia el títol de la finestra

        #self.showMaximized()   #Maximitza la finestra
        self.setMinimumSize(250, 400)   #Mida mínima (x, y)
        self.setMaximumSize(250, 400)   #Mida màxima de la finestra



app = QApplication(sys.argv)

_menu_prin = Menu_prin()

_menu_prin.show()


app.exec_()
