
import os
#from threading import Thread
from rtlsdr import RtlSdr

import time
import numpy as np 

#from PyQt5.QtGui import QFont
from PyQt5 import QtWidgets,QtCore
#from PySide2.QtUiTools import QUiLoader
#import pyqtgraph as pg

import sys

ui_cwd = os.getcwd()

sys.path.append( str(ui_cwd)+"/ui/" )
from SDR import Ui_MainWindow
from RTL_sdr import Rtlsdr

today = time.strftime("%Y-%m-%d", time.localtime())
file_name = "./log/"+str(today) + "_sdr_log.txt"

class SDR_gui(QtWidgets.QMainWindow):
    def __init__(self):
        self.rtlsdr = None
        self.currentTime = None
        self.signal_data = np.array([])
        super().__init__()
        self.ui = Ui_MainWindow()
        # 初始化界面

        self.ui.setupUi(self)

        self.ui.signal.setLimits(xMin=24,xMax=1700,yMin=-90,yMax=90)
        self.ui.signal.setXRange(min=24,max=1700)
        self.ui.signal.setYRange(min=-90,max=90)
        self.ui.signal.setBackground('w')
        self.ui.signal.setLabel("bottom","Frequency (MHz)",size="10pt")
        self.ui.signal.setLabel("left","Relative Power (dB)",size="10pt")
        self.ui.signal.showGrid(x=True, y=True)

        self.ui.start.clicked.connect(self.start_sdr)
        self.ui.start.clicked.connect(self.obs_signal)
        self.ui.dis_start.clicked.connect(self.end_sdr)
        self.ui.clear.clicked.connect(self.clear_data)
        self.ui.save.clicked.connect(self.save_file)

        with open(file_name,"a") as f:
            f.write("####################\n\n####################\n")
            f.write("SDR GUI start\n")
            f.write("####################\n\n####################\n")

    def start_sdr(self):
        if self.rtlsdr == None:
            self.rtlsdr = Rtlsdr()
        else:
            pass
        
        self.rtlsdr.sdr.gain = float(self.ui.gain.toPlainText())
        self.rtlsdr.sdr.center_freq = float( self.ui.cfreq.toPlainText() )*1e6
        self.obstime = int( np.round(float( self.ui.obstime.toPlainText() )*60) )
        temp_fftsize = int( self.ui.fftsize.toPlainText() )

        if temp_fftsize == 0:
            self.fftsize = 2
        else:
            self.fftsize = int( 2**np.round( np.log2(temp_fftsize) ) )


        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("##########\n")
            f.write(f"{temp_time}: SDR start work\n")
            f.write("##########\n")

        if self.ui.sample_rate.currentText() == '0.24Mhz':
            self.rtlsdr.sdr.sample_rate = 0.24e6
        elif self.ui.sample_rate.currentText() == '0.3Mhz':
            self.rtlsdr.sdr.sample_rate = 0.3e6
        elif self.ui.sample_rate.currentText() == '0.96Mhz':
            self.rtlsdr.sdr.sample_rate = 0.96e6
        elif self.ui.sample_rate.currentText() == '1.2Mhz':
            self.rtlsdr.sdr.sample_rate = 1.2e6
        elif self.ui.sample_rate.currentText() == '1.44Mhz':
            self.rtlsdr.sdr.sample_rate = 1.44e6
        elif self.ui.sample_rate.currentText() == '1.6Mhz':
            self.rtlsdr.sdr.sample_rate = 1.6e6
        elif self.ui.sample_rate.currentText() == '1.8Mhz':
            self.rtlsdr.sdr.sample_rate = 1.8e6
        elif self.ui.sample_rate.currentText() == '2.4Mhz':
            self.rtlsdr.sdr.sample_rate = 2.4e6
        else:
            self.rtlsdr.sdr.sample_rate = 3.2e6
        
        with open(file_name,"a") as f:
                f.write("##########\n")
                f.write(f"Center frequency:{self.rtlsdr.sdr.center_freq} Hz.\n")
                f.write(f"Sample rate: {self.rtlsdr.sdr.sample_rate} Hz.\n")
                f.write(f"Gain: {self.rtlsdr.sdr.gain} dB.\n")
                f.write(f"FFTsize: {self.fftsize} \n")
                f.write(f"Observation time: {self.obstime}.\n")
                f.write("##########\n")

        print(f"sdr cfreq is {self.rtlsdr.sdr.center_freq}, sdr fftsize is {self.fftsize}, sdr gain is {self.rtlsdr.sdr.gain}, sdr sample rate is {self.rtlsdr.sdr.sample_rate}.")
        print(f"obs time is {self.obstime}.")


    def plot_data(self):

        time.sleep(0.2)
        
        freql,datal = self.rtlsdr.obs_1s(self.rtlsdr.sdr.center_freq,self.rtlsdr.sdr.sample_rate,self.fftsize,self.rtlsdr.sdr.gain)
        if len(self.signal_data) == 0:
            self.signal_data = datal
        else:
            self.signal_data = np.hstack( (self.signal_data,datal) )
        xl = freql/1e6
        yl = 10*np.log10(datal)

        self.ui.signal.clear()
        xmin = np.min(xl) - 0.05*self.rtlsdr.sdr.sample_rate/1e6
        xmax = np.max(xl) + 0.05*self.rtlsdr.sdr.sample_rate/1e6
        ymin = np.min(yl) - 0.5
        ymax =np.max(yl) + 0.5
        self.ui.signal.setLimits(xMin=xmin,xMax=xmax,yMin=ymin,yMax=ymax)
        self.ui.signal.setXRange(min=xmin,max=xmax)
        self.ui.signal.setYRange(min=ymin,max=ymax)

        self.ui.signal.plot().setData( xl,yl,pen=(21,101,192) )
        
    def obs_signal(self):
        self.currentTime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        self.timer = QtCore.QTimer(self)
        self.timer2 = QtCore.QTimer(self)
        self.timer.timeout.connect(self.plot_data)
        self.ui.stop.clicked.connect(self.broken_timer)
        self.timer.start(1050)
        self.timer2.singleShot( int(250+1050*(self.obstime*2-2)),self.timer.stop )

    def broken_timer(self):
        self.timer.stop()
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write(f"SDR stopped at {temp_time}.\n##########\n")


    def clear_data(self):
        self.currentTime = None
        self.signal_data = np.array([])
        self.ui.signal.clear()
        print("data cleaned.")
        with open(file_name,"a") as f:
            f.write("data cleaned.\n")

    def save_file(self):
        file_path = self.ui.datapath.toPlainText()
        file_name_data = str(self.rtlsdr.sdr.center_freq/1e6) + "MHz_" + str(self.rtlsdr.sdr.sample_rate/1e6) + "MHz_" + str(self.fftsize) + "fftsize_" + str(self.obstime) + "s_" + str(self.rtlsdr.sdr.gain) + "dB.txt"
        print(file_path+file_name_data)
        if self.currentTime == None:
            print("No data, pass...\n")
        else:
            self.signal_data = self.signal_data.reshape(self.fftsize,int(len(self.signal_data)/self.fftsize))
            print("data saved!")
            with open(file_name,"a") as f:
                f.write("data saved.\n")
            np.savetxt( str(file_path+file_name_data),self.signal_data )
            

    def end_sdr(self):
        print("SDR Disconnected.")
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("####################\n")
            f.write(f"{temp_time}: SDR disconnect.\n")
            f.write("####################\n")
        self.rtlsdr.sdr.close()

if __name__ == '__main__':
    app = QtWidgets.QApplication([sys.argv])
    mainw = SDR_gui()
    mainw.show()
    sys.exit(app.exec_())



