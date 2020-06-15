# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'plot_window_ui.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindowPlot(object):
    def setupUi(self, MainWindowPlot):
        MainWindowPlot.setObjectName("MainWindowPlot")
        MainWindowPlot.resize(800, 600)
        MainWindowPlot.setIconSize(QtCore.QSize(24, 24))
        self.centralwidget = QtWidgets.QWidget(MainWindowPlot)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        MainWindowPlot.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindowPlot)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        MainWindowPlot.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindowPlot)
        self.statusbar.setEnabled(True)
        self.statusbar.setObjectName("statusbar")
        MainWindowPlot.setStatusBar(self.statusbar)
        self.toolBar = QToolBarExpanding(MainWindowPlot)
        self.toolBar.setIconSize(QtCore.QSize(20, 20))
        self.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        self.toolBar.setObjectName("toolBar")
        MainWindowPlot.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionPan = QtWidgets.QAction(MainWindowPlot)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/plot/pan"), QtGui.QIcon.Normal, QtGui.QIcon.On)
        self.actionPan.setIcon(icon)
        self.actionPan.setObjectName("actionPan")
        self.actionZoom = QtWidgets.QAction(MainWindowPlot)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/plot/zoom"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionZoom.setIcon(icon1)
        self.actionZoom.setObjectName("actionZoom")
        self.actionConnectors = QtWidgets.QAction(MainWindowPlot)
        self.actionConnectors.setCheckable(True)
        self.actionConnectors.setChecked(False)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/connectors"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionConnectors.setIcon(icon2)
        self.actionConnectors.setObjectName("actionConnectors")
        self.actionCoords = QtWidgets.QAction(MainWindowPlot)
        self.actionCoords.setCheckable(True)
        self.actionCoords.setChecked(True)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/plot/point"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionCoords.setIcon(icon3)
        self.actionCoords.setObjectName("actionCoords")
        self.actionAuto = QtWidgets.QAction(MainWindowPlot)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/plot/autozoom"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionAuto.setIcon(icon4)
        self.actionAuto.setObjectName("actionAuto")
        self.actionReplot = QtWidgets.QAction(MainWindowPlot)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/plot/refresh_plot"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionReplot.setIcon(icon5)
        self.actionReplot.setObjectName("actionReplot")
        self.actionRuler = QtWidgets.QAction(MainWindowPlot)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap(":/plot/ruler"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionRuler.setIcon(icon6)
        self.actionRuler.setObjectName("actionRuler")
        self.toolBar.addAction(self.actionPan)
        self.toolBar.addAction(self.actionAuto)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionCoords)
        self.toolBar.addAction(self.actionConnectors)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionRuler)

        self.retranslateUi(MainWindowPlot)
        self.actionAuto.triggered.connect(MainWindowPlot.auto_scale)
        self.actionConnectors.triggered['bool'].connect(MainWindowPlot.set_show_connectors)
        self.actionCoords.triggered['bool'].connect(MainWindowPlot.set_position_track)
        self.actionPan.triggered.connect(MainWindowPlot.pan)
        self.actionZoom.triggered.connect(MainWindowPlot.zoom)
        self.actionReplot.triggered.connect(MainWindowPlot.replot)
        QtCore.QMetaObject.connectSlotsByName(MainWindowPlot)

    def retranslateUi(self, MainWindowPlot):
        _translate = QtCore.QCoreApplication.translate
        MainWindowPlot.setWindowTitle(_translate("MainWindowPlot", "MainWindow"))
        self.toolBar.setWindowTitle(_translate("MainWindowPlot", "toolBar"))
        self.actionPan.setText(_translate("MainWindowPlot", "Help"))
        self.actionPan.setShortcut(_translate("MainWindowPlot", "P"))
        self.actionZoom.setText(_translate("MainWindowPlot", "Zoom"))
        self.actionZoom.setToolTip(_translate("MainWindowPlot", "Zoom control"))
        self.actionZoom.setShortcut(_translate("MainWindowPlot", "Z"))
        self.actionConnectors.setText(_translate("MainWindowPlot", "Pins"))
        self.actionConnectors.setToolTip(_translate("MainWindowPlot", "Show connectors pins for selected qcomponents"))
        self.actionConnectors.setShortcut(_translate("MainWindowPlot", "C"))
        self.actionCoords.setText(_translate("MainWindowPlot", "Get point"))
        self.actionCoords.setToolTip(_translate("MainWindowPlot", "Click for position --- Enable this to click on the plot and log the (x,y) position"))
        self.actionCoords.setShortcut(_translate("MainWindowPlot", "P"))
        self.actionAuto.setText(_translate("MainWindowPlot", "Autoscale"))
        self.actionAuto.setToolTip(_translate("MainWindowPlot", "Auto Zoom"))
        self.actionAuto.setShortcut(_translate("MainWindowPlot", "A"))
        self.actionReplot.setText(_translate("MainWindowPlot", "Replot"))
        self.actionReplot.setShortcut(_translate("MainWindowPlot", "Ctrl+R"))
        self.actionRuler.setText(_translate("MainWindowPlot", "Ruler"))
        self.actionRuler.setToolTip(_translate("MainWindowPlot", "Activate the ruler"))
from .widgets.bases.expanding_toolbar import QToolBarExpanding
from . import main_window_rc_rc
