# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'location.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Set_location(object):
    def setupUi(self, Set_location):
        Set_location.setObjectName("Set_location")
        Set_location.resize(720, 540)
        Set_location.setMinimumSize(QtCore.QSize(720, 540))
        Set_location.setMaximumSize(QtCore.QSize(720, 540))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        Set_location.setFont(font)
        self.input_location = QtWidgets.QTextEdit(Set_location)
        self.input_location.setGeometry(QtCore.QRect(180, 30, 400, 40))
        self.input_location.setObjectName("input_location")
        self.label_altitude = QtWidgets.QLabel(Set_location)
        self.label_altitude.setGeometry(QtCore.QRect(60, 90, 100, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_altitude.setFont(font)
        self.label_altitude.setAlignment(QtCore.Qt.AlignCenter)
        self.label_altitude.setObjectName("label_altitude")
        self.label_longitude = QtWidgets.QLabel(Set_location)
        self.label_longitude.setGeometry(QtCore.QRect(60, 170, 100, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_longitude.setFont(font)
        self.label_longitude.setAlignment(QtCore.Qt.AlignCenter)
        self.label_longitude.setObjectName("label_longitude")
        self.label_height = QtWidgets.QLabel(Set_location)
        self.label_height.setGeometry(QtCore.QRect(60, 250, 100, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_height.setFont(font)
        self.label_height.setAlignment(QtCore.Qt.AlignCenter)
        self.label_height.setObjectName("label_height")
        self.input_altitude = QtWidgets.QTextEdit(Set_location)
        self.input_altitude.setGeometry(QtCore.QRect(180, 90, 400, 40))
        self.input_altitude.setObjectName("input_altitude")
        self.input_longitude = QtWidgets.QTextEdit(Set_location)
        self.input_longitude.setGeometry(QtCore.QRect(180, 170, 400, 40))
        self.input_longitude.setObjectName("input_longitude")
        self.input_height = QtWidgets.QTextEdit(Set_location)
        self.input_height.setGeometry(QtCore.QRect(180, 250, 400, 40))
        self.input_height.setObjectName("input_height")
        self.label_timeoffset = QtWidgets.QLabel(Set_location)
        self.label_timeoffset.setGeometry(QtCore.QRect(60, 330, 121, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_timeoffset.setFont(font)
        self.label_timeoffset.setAlignment(QtCore.Qt.AlignCenter)
        self.label_timeoffset.setObjectName("label_timeoffset")
        self.select_timezone = QtWidgets.QComboBox(Set_location)
        self.select_timezone.setGeometry(QtCore.QRect(190, 330, 100, 40))
        self.select_timezone.setObjectName("select_timezone")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.select_timezone.addItem("")
        self.buttonBox = QtWidgets.QDialogButtonBox(Set_location)
        self.buttonBox.setGeometry(QtCore.QRect(500, 440, 171, 81))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.label_location = QtWidgets.QLabel(Set_location)
        self.label_location.setGeometry(QtCore.QRect(60, 30, 100, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_location.setFont(font)
        self.label_location.setAlignment(QtCore.Qt.AlignCenter)
        self.label_location.setObjectName("label_location")
        self.label_rotator = QtWidgets.QLabel(Set_location)
        self.label_rotator.setGeometry(QtCore.QRect(300, 330, 121, 40))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_rotator.setFont(font)
        self.label_rotator.setAlignment(QtCore.Qt.AlignCenter)
        self.label_rotator.setObjectName("label_rotator")
        self.input_port = QtWidgets.QTextEdit(Set_location)
        self.input_port.setGeometry(QtCore.QRect(430, 330, 211, 40))
        self.input_port.setObjectName("input_port")

        self.retranslateUi(Set_location)
        QtCore.QMetaObject.connectSlotsByName(Set_location)

    def retranslateUi(self, Set_location):
        _translate = QtCore.QCoreApplication.translate
        Set_location.setWindowTitle(_translate("Set_location", "Location Setting"))
        self.input_location.setHtml(_translate("Set_location", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600; font-style:italic;\">\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Boston</p></body></html>"))
        self.label_altitude.setText(_translate("Set_location", "Altitude"))
        self.label_longitude.setText(_translate("Set_location", "Longitude"))
        self.label_height.setText(_translate("Set_location", "Height"))
        self.input_altitude.setHtml(_translate("Set_location", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600; font-style:italic;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">42.36</p></body></html>"))
        self.input_longitude.setHtml(_translate("Set_location", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600; font-style:italic;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">-71.06</p></body></html>"))
        self.input_height.setHtml(_translate("Set_location", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600; font-style:italic;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">43</p></body></html>"))
        self.label_timeoffset.setText(_translate("Set_location", "Time Offset"))
        self.select_timezone.setCurrentText(_translate("Set_location", "-12"))
        self.select_timezone.setItemText(0, _translate("Set_location", "-12"))
        self.select_timezone.setItemText(1, _translate("Set_location", "-11"))
        self.select_timezone.setItemText(2, _translate("Set_location", "-10"))
        self.select_timezone.setItemText(3, _translate("Set_location", "-9"))
        self.select_timezone.setItemText(4, _translate("Set_location", "-8"))
        self.select_timezone.setItemText(5, _translate("Set_location", "-7"))
        self.select_timezone.setItemText(6, _translate("Set_location", "-6"))
        self.select_timezone.setItemText(7, _translate("Set_location", "-5"))
        self.select_timezone.setItemText(8, _translate("Set_location", "-4"))
        self.select_timezone.setItemText(9, _translate("Set_location", "-3"))
        self.select_timezone.setItemText(10, _translate("Set_location", "-2"))
        self.select_timezone.setItemText(11, _translate("Set_location", "-1"))
        self.select_timezone.setItemText(12, _translate("Set_location", "0"))
        self.select_timezone.setItemText(13, _translate("Set_location", "1"))
        self.select_timezone.setItemText(14, _translate("Set_location", "2"))
        self.select_timezone.setItemText(15, _translate("Set_location", "3"))
        self.select_timezone.setItemText(16, _translate("Set_location", "4"))
        self.select_timezone.setItemText(17, _translate("Set_location", "5"))
        self.select_timezone.setItemText(18, _translate("Set_location", "6"))
        self.select_timezone.setItemText(19, _translate("Set_location", "7"))
        self.select_timezone.setItemText(20, _translate("Set_location", "8"))
        self.select_timezone.setItemText(21, _translate("Set_location", "9"))
        self.select_timezone.setItemText(22, _translate("Set_location", "10"))
        self.select_timezone.setItemText(23, _translate("Set_location", "11"))
        self.select_timezone.setItemText(24, _translate("Set_location", "12"))
        self.label_location.setText(_translate("Set_location", "Location"))
        self.label_rotator.setText(_translate("Set_location", "Rotator"))
        self.input_port.setHtml(_translate("Set_location", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:600; font-style:italic;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">/dev/ttyUSB0</p></body></html>"))
