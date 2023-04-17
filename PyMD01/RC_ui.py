
import os
#from multiprocessing import Process
from threading import Thread

import astropy.units as u
from  astropy.time import  Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon

import time
import numpy as np 
import matplotlib.pyplot as plt 

from astropy.visualization import astropy_mpl_style, quantity_support 
plt.style.use(astropy_mpl_style) 
quantity_support()

from PyQt5.QtGui import QFont
from PyQt5 import QtWidgets,QtCore
#from PySide2.QtUiTools import QUiLoader
import pyqtgraph as pg

import sys

ui_cwd = os.getcwd()

sys.path.append( str(ui_cwd)+"/ui/" )
from PyMd01_RC import Ui_MainWindow
from PyMd01 import RotatorController
from location import Ui_Set_location
#from SDR_ui import SDR_gui

today = time.strftime("%Y-%m-%d", time.localtime())
file_name = "./log/"+str(today) + "_rotator_log.txt"

Boston = EarthLocation(lat=42.36*u.deg, lon=-71.06*u.deg, height=43*u.m) 
timezone = -5

az_origin = np.arange(180,-210,-30)
az_show = ['N','330','300','W','240','210','S','150','120','E','60','30','N']


def calculate_sky():
    """
    Computing the Galactic plane position and the positions of the Sun and Moon
    """
    xl = []
    yl = []
    ctime1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    ctime = Time(ctime1)-timezone*u.hour  # Eastern Daylight Time 

    lrange = np.arange(0,360,10)
    brange = np.linspace(0,0,36)
    loc_sun = get_sun(ctime).transform_to( AltAz(obstime=ctime,location=Boston) )
    loc_moon = get_moon(ctime).transform_to( AltAz(obstime=ctime,location=Boston) )
    sun_az = loc_sun.az.value - 180
    sun_alt = loc_sun.alt.value

    moon_az = loc_moon.az.value - 180
    moon_alt = loc_moon.alt.value


    prange = SkyCoord(lrange*u.degree, brange*u.degree, frame='galactic') 
    praltaz = prange.transform_to(AltAz(obstime=ctime,location=Boston))
    for i in range(len(praltaz.alt.value)):
        temp_az = praltaz.az.value[i] - 180
        xl.append(temp_az)
        yl.append(praltaz.alt.value[i]) 
    return ctime1,xl,yl,lrange,brange,sun_az,sun_alt,moon_az,moon_alt


def set_calculate_sky(Alt,Lon,Hei,name,tiz):
    """
    Calculating the position of the Galactic plane and the Sun and Moon.
    """
    name = EarthLocation(lat=Alt*u.deg, lon=Lon*u.deg, height=Hei*u.m) 
    xl = []
    yl = []
    ctime1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    ctime = Time(ctime1)-tiz*u.hour  # Eastern Daylight Time 

    lrange = np.arange(0,360,10)
    brange = np.linspace(0,0,36)
    loc_sun = get_sun(ctime).transform_to( AltAz(obstime=ctime,location=name) )
    loc_moon = get_moon(ctime).transform_to( AltAz(obstime=ctime,location=name) )
    sun_az = loc_sun.az.value - 180
    sun_alt = loc_sun.alt.value

    moon_az = loc_moon.az.value - 180
    moon_alt = loc_moon.alt.value


    prange = SkyCoord(lrange*u.degree, brange*u.degree, frame='galactic') 
    praltaz = prange.transform_to(AltAz(obstime=ctime,location=name))
    for i in range(len(praltaz.alt.value)):
        temp_az = praltaz.az.value[i] - 180
        xl.append(temp_az)
        yl.append(praltaz.alt.value[i]) 
    return ctime1,xl,yl,lrange,brange,sun_az,sun_alt,moon_az,moon_alt


def AZEL2ICRS(AZ,EL,Alt,Lon,Hei,tiz,name):
    name = EarthLocation(lat=Alt*u.deg, lon=Lon*u.deg, height=Hei*u.m) 
    ctime1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    ctime = Time(ctime1)-tiz*u.hour  # Eastern Daylight Time 
    obj_coord = SkyCoord( AltAz(az = AZ*u.degree,
    alt = EL*u.degree,
    obstime=ctime,location=name) )
    obj_rd = obj_coord.transform_to('icrs')
    return obj_rd


def Gal2ICRS(L,B):
    obj_coord = SkyCoord(float(L)*u.degree, float(B)*u.degree, frame='galactic')
    obj_rd = obj_coord.icrs

    return obj_rd


################################################################################
################################################################################

class TrackDialog(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Notice!")

        self.setFixedSize(480,270)
        self.layout = QtWidgets.QVBoxLayout()
        message = QtWidgets.QLabel("Input AZALT coordiantes, tracking is useless.\n\nOnce you start tracking an object, you'd better\n not change anything untill tracking stop.")
        message.setFont(QFont("Times", 14, QFont.Bold))
        message.setAlignment(QtCore.Qt.AlignCenter)
        self.layout.addWidget(message)
        #self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)


class setting_location(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Set_location()
        self.ui.setupUi(self)
        self.ui.select_timezone.setCurrentIndex(7)
        self.ui.buttonBox.accepted.connect(self.get_info)
        self.ui.buttonBox.accepted.connect(self.close)
        self.ui.buttonBox.rejected.connect(self.close)

    def get_info(self):
        self.loc = str( self.ui.input_location.toPlainText() )
        self.alt = float( self.ui.input_altitude.toPlainText() )
        self.lon = float( self.ui.input_longitude.toPlainText() )
        self.hei = float( self.ui.input_height.toPlainText() )
        self.tiz = float( self.ui.select_timezone.currentText() )
        self.port = str( self.ui.input_port.toPlainText() )

        with open(file_name,"a") as f:
            f.write("####################\n\n####################\nStart Rotator Controller\n")
            f.write("Location information:\n")
            f.write(f"Location:{self.loc}\n")
            f.write(f"Latitude(deg):{self.alt}\n")
            f.write(f"Longitude(deg):{self.lon}\n")
            f.write(f"Altitude(m):{self.hei}\n")
            f.write(f"Timezone(h):{self.tiz}\n")
            f.write(f"Rotator Port:{self.port}\n")
            f.write("##########\n")
        print(self.loc,self.alt,self.lon,self.hei,self.tiz,self.port)
        

class Rc_Window(QtWidgets.QMainWindow):

    def __init__(self,loc_name,loc_alt,loc_lon,loc_hei,loc_tiz,port):
        super().__init__()
        # Importing the user interface using a UI file.
        self.ui = Ui_MainWindow()
        # initial

        self.ui.setupUi(self)
        
        self.file_write_count = 0
        self.tra = None
        self.tdec = None

        self.loc_name = loc_name
        self.loc_alt = loc_alt
        self.loc_lon = loc_lon
        self.loc_hei = loc_hei
        self.loc_tiz = loc_tiz
        self.port = port


        self.point_az = 0
        self.point_alt = 0

        # initial sky
        self.ui.sky_area.setLimits(xMin=-190,xMax=190,yMin=0,yMax=95)
        self.ui.sky_area.setXRange(min=-190,max=190)
        self.ui.sky_area.setYRange(min=0,max=91)
        self.ui.sky_area.setBackground('w')
        self.ui.sky_area.setLabel("bottom","Azimuth (deg)",size="10pt")
        self.ui.sky_area.setLabel("left","Altitude (deg)",size="10pt")
        self.ui.sky_area.showGrid(x=True, y=True)

        az_xticks = [ [i,j] for i,j in zip(az_origin,az_show) ]
        self.ui.sky_area.getAxis('bottom').setTicks([az_xticks])
        
        # Mutil process
        sky_process = Thread(target = self.ui.start_sky.clicked.connect(self.timer_start))
        rc_process = Thread(target = self.ui.rc_start.clicked.connect(self.start_rotator))
        select_obj = Thread( target=self.selection() )
        track_process = Thread( target=self.ui.tracking.clicked.connect(self.tracking_timer) )
        sky_process.start()
        rc_process.start()
        select_obj.start()
        track_process.start()
        
        # init buttons
        self.ui.stop_sky.clicked.connect(self.timer_end)
        
        self.ui.rc_stop.clicked.connect(self.stop_rotator)
        
        # init text input
        self.ui.coord1_input.setText("Input coordinate 1")

        self.ui.coord2_input.setText("Input coordinate 2")

        self.ui.track_time.setText("Example: tracking time in mins 2.1")

        #self.ui.tra
        track_process.join()
        select_obj.join()
        sky_process.join()
        rc_process.join()

############################################################

    def tracking_radec(self):
        if self.ui.obj_name.currentText() == "Sun":
            t_hours = float( self.ui.track_time.toPlainText() )
            self.rc.tracking_sun(hours=t_hours)
            temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            with open(file_name,"a") as f:
                f.write("Tracking target:Sun\n")
                f.write(f"Start time:{temp_time}\n")
                f.write(f"Tracking hours:{t_hours}\n")
        elif self.ui.obj_name.currentText() == "Moon":
            t_hours = float( self.ui.track_time.toPlainText() )
            self.rc.tracking_moon(hours=t_hours)
            temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            with open(file_name,"a") as f:
                f.write("Tracking target:Moon\n")
                f.write(f"Start time:{temp_time}\n")
                f.write(f"Tracking hours:{t_hours}\n")
        else:
            if self.tra == None:
                print("\nNo tracking coordinates.Plz Goto target first.\n")
            else:
                t_hours = float( self.ui.track_time.toPlainText() )
                self.rc.tracking( self.tra,self.tdec,hours=t_hours )
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write(f"Tracking target:Ra={self.tra};Dec={self.tdec}\n")
                    f.write(f"Start time:{temp_time}\n")
                    f.write(f"Tracking hours:{t_hours}\n")

    def selection(self):
        self.ui.obj_name.currentTextChanged.connect(self.name_obj)
        self.ui.select_coords.currentTextChanged.connect(self.name_coords)
        self.ui.goto_2.clicked.connect( self.goto_obj )

    def goto_obj(self):
        if self.ui.obj_name.currentText() == "NA":
            if self.ui.select_coords.currentText() == "AZEL - AZ_ALT":
                print("Input AZEL coordinates.")
                input_az = self.ui.coord1_input.toPlainText()
                input_el = self.ui.coord2_input.toPlainText()
                self.rc.AZEL_goto_100( float(input_az),float(input_el) )
                tracking_rd = AZEL2ICRS( float(input_az),float(input_el),self.loc_alt,self.loc_lon,self.loc_hei,self.loc_tiz,self.loc_name )
                self.point_az = float(input_az)
                self.point_alt = float(input_el)
                self.tra = str(tracking_rd.ra)
                self.tdec = str(tracking_rd.dec)
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write("**********\n")
                    f.write(f"Time:{temp_time}. ")
                    f.write("AZEL\n")
                    f.write(f"Target goto Az={self.point_az},Alt={self.point_alt}\n")
                    f.write(f"Target goto Ra={self.tra},Dec={self.tdec}\n")
                    f.write("**********\n")

            elif self.ui.select_coords.currentText() == "ICRS - RA_DEC":
                print("Input ICRS coordinates.")
                input_ra = self.ui.coord1_input.toPlainText()
                input_dec = self.ui.coord2_input.toPlainText()
                self.point_az,self.point_alt = self.rc.RaDec_goto(str(input_ra),str(input_dec))
                self.tra = str(input_ra)
                self.tdec = str(input_dec)
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write("**********\n")
                    f.write(f"Time:{temp_time}. ")
                    f.write("ICRS\n")
                    f.write(f"Target goto Az={self.point_az},Alt={self.point_alt}\n")
                    f.write(f"Target goto Ra={self.tra},Dec={self.tdec}\n")
                    f.write("**********\n")

            elif self.ui.select_coords.currentText() == "Galatic - l_b":
                print("Input Galatic coordinates.")
                input_l = self.ui.coord1_input.toPlainText()
                input_b = self.ui.coord2_input.toPlainText()
                tracking_rd = Gal2ICRS(input_l,input_b)
                self.point_az,self.point_alt = self.rc.RaDec_goto(str(tracking_rd.ra),str(tracking_rd.dec))
                self.tra = str(tracking_rd.ra)
                self.tdec = str(tracking_rd.dec)
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write("**********\n")
                    f.write(f"Time:{temp_time}. ")
                    f.write("Galatic\n")
                    f.write(f"Target goto l={input_l},b={input_b}\n")
                    f.write(f"Target goto Ra={self.tra},Dec={self.tdec}\n")
                    f.write("**********\n")

            else:
                print("......")
        else:
            if self.ui.obj_name.currentText() == "Sun":
                self.rc.goto_sun()
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write("**********\n")
                    f.write(f"Time:{temp_time}. ")
                    f.write("Target goto Sun\n")
                    f.write("**********\n")
            elif self.ui.obj_name.currentText() == "Moon":
                self.rc.goto_moon()
                temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                with open(file_name,"a") as f:
                    f.write("**********\n")
                    f.write(f"Time:{temp_time}. ")
                    f.write("Target goto Moon\n")
                    f.write("**********\n")
            elif self.ui.obj_name.currentText() == "M31":
                print("Goto M31.....")
            else:
                print("GoTo by coords...")


    def name_obj(self):
        if self.ui.obj_name.currentText() == 'Sun':
            print("You selected Sun.\n")
            #self.set_move_button_dis()
            self.ui.coord1_input.setDisabled(True)
            self.ui.coord2_input.setDisabled(True)

        elif self.ui.obj_name.currentText() == 'Moon':
            print("You selected Moon.\n")
            #self.set_move_button_dis()
            self.ui.coord1_input.setDisabled(True)
            self.ui.coord2_input.setDisabled(True)

        elif self.ui.obj_name.currentText() == 'M31':
            print("You selected M31.\n")
            #self.set_move_button_dis()
            self.ui.coord1_input.setDisabled(True)
            self.ui.coord2_input.setDisabled(True)

        elif self.ui.obj_name.currentText() == 'NA':
            #self.set_move_button_en()
            self.ui.coord1_input.setDisabled(False)
            self.ui.coord2_input.setDisabled(False)
            print( "No selected object.\n" )

        else:
            #self.set_move_button_en()
            self.ui.coord1_input.setDisabled(True)
            self.ui.coord2_input.setDisabled(True)


    def name_coords(self):
        if self.ui.select_coords.currentText() == "AZEL - AZ_ALT":

            self.ui.coord1_input.setDisabled(False)
            self.ui.coord2_input.setDisabled(False)
            self.ui.coord1_input.setReadOnly(False)
            self.ui.coord2_input.setReadOnly(False)
            self.ui.coord1_input.setText("Example: 220.34")
            self.ui.coord2_input.setText("Example: 30.34")
            print("Current Coordinates system is AZEL.")

        elif self.ui.select_coords.currentText() == "ICRS - RA_DEC":

            self.ui.coord1_input.setDisabled(False)
            self.ui.coord2_input.setDisabled(False)
            self.ui.coord1_input.setReadOnly(False)
            self.ui.coord2_input.setReadOnly(False)
            self.ui.coord1_input.setText("Example: 02h30m22s")
            self.ui.coord2_input.setText("Example: 20d45m11s")
            print("Current Coordinates system is ICRS.")
        elif self.ui.select_coords.currentText() == "Galatic - l_b":

            self.ui.coord1_input.setDisabled(False)
            self.ui.coord2_input.setDisabled(False)
            self.ui.coord1_input.setReadOnly(False)
            self.ui.coord2_input.setReadOnly(False)
            self.ui.coord1_input.setText("Example: 230.234")
            self.ui.coord2_input.setText("Example: 3.75")
            print("Current Coordinates system is Galatic.")
        else:
            
            self.ui.coord1_input.setDisabled(False)
            self.ui.coord2_input.setDisabled(False)
            self.ui.coord1_input.setReadOnly(True)
            self.ui.coord2_input.setReadOnly(True)
            self.ui.coord1_input.setText("No coordinates system selected.")
            self.ui.coord2_input.setText("No coordinates system selected.")
            print("You not select a coordinates system. Using object name...\n")

    def start_rotator(self):
        
        self.rc = RotatorController(self.port,loclon=self.loc_lon,loclat=self.loc_alt,locheight=self.loc_hei,utcoffset=self.loc_tiz)
        if self.rc.ser.is_open:
            self.rc.get_loc_100()
            self.ui.label_rc.setText('Rotator Connected.')
            print( "The rotator has connected.\n" )
            print( f"Current pointing location is: az = {self.rc.az_100}; el = {self.rc.el_100}.\n" )
            
            dlg = TrackDialog()
            dlg.exec_()
        else:
            print( "Need connect the rotator.\n" )
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f"Time:{temp_time}. ")
            f.write("Rotator connected\n")
            f.write("**********\n")

        self.ui.up.clicked.connect( self.rc.move_up )
        self.ui.left_up.clicked.connect( self.rc.move_left_up )
        self.ui.left.clicked.connect( self.rc.move_left )
        self.ui.left_down.clicked.connect( self.rc.move_left_down )
        self.ui.down.clicked.connect( self.rc.move_down )
        self.ui.right_down.clicked.connect( self.rc.move_right_down )
        self.ui.right.clicked.connect( self.rc.move_right )
        self.ui.right_up.clicked.connect( self.rc.move_right_up )

        self.ui.stop.clicked.connect( self.rc.stop )
        self.ui.stop.clicked.connect( self.wirte_stop_log )
        self.ui.stop.clicked.connect( self.rc.get_loc_100 )
        self.ui.stop.clicked.connect( self.update_point )
        self.ui.stop.clicked.connect( self.tracking_timer_end )

    def update_point(self):
        self.point_az = self.rc.az_100
        self.point_alt = self.rc.el_100

    def wirte_stop_log(self):
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f" Rotator stopped by usr at {temp_time}. Az = {self.point_az}; Alt = {self.point_alt}.\n")
            f.write("**********\n")

    def stop_rotator(self):
        self.rc.close_comm()
        print( "The rotator is closing...\n" )
        self.ui.label_rc.setText('Rotator Disconnected.')
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time.sleep(0.1)
        self.point_az = 0
        self.point_alt = 0
        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f"Time:{temp_time}.\n")
            f.write("Rotator Disconnected.\n")
            f.write("**********\n")

    def tracking_obj(self):
        if self.ui.obj_name.currentText() == "NA":
            if self.ui.select_coords.currentText() == "AZEL - AZ_ALT":
                input_az = self.ui.coord1_input.toPlainText()
                input_el = self.ui.coord2_input.toPlainText()
                tracking_rd = AZEL2ICRS( float(input_az),float(input_el),self.loc_alt,self.loc_lon,self.loc_hei,self.loc_tiz,self.loc_name )
                self.rc.RaDec_goto( str(tracking_rd.ra),str(tracking_rd.dec) )
                self.point_az = float(input_az)
                self.point_alt = float(input_el)
                self.tra = str(tracking_rd.ra)
                self.tdec = str(tracking_rd.dec)
                if self.file_write_count == 0:
                    with open(file_name,"a") as f:
                        f.write(f"Tracking target(AZEL): Az = {self.point_az}; Alt = {self.point_alt}. Ra = {self.tra}; Dec = {self.tdec}.\n ")
                        self.file_write_count += 1
                else:
                    pass

            elif self.ui.select_coords.currentText() == "ICRS - RA_DEC":
                input_ra = self.ui.coord1_input.toPlainText()
                input_dec = self.ui.coord2_input.toPlainText()
                self.point_az,self.point_alt = self.rc.RaDec_goto(str(input_ra),str(input_dec))
                self.tra = str(input_ra)
                self.tdec = str(input_dec)

                if self.file_write_count == 0:
                    with open(file_name,"a") as f:
                        f.write(f"Tracking target(ICRS): Az = {self.point_az}; Alt = {self.point_alt}. Ra = {self.tra}; Dec = {self.tdec}.\n ")
                        self.file_write_count += 1

            elif self.ui.select_coords.currentText() == "Galatic - l_b":
                input_l = self.ui.coord1_input.toPlainText()
                input_b = self.ui.coord2_input.toPlainText()
                tracking_rd = Gal2ICRS(input_l,input_b)
                self.point_az,self.point_alt = self.rc.RaDec_goto(str(tracking_rd.ra),str(tracking_rd.dec))
                self.tra = str(tracking_rd.ra)
                self.tdec = str(tracking_rd.dec)
                if self.file_write_count == 0:
                    with open(file_name,"a") as f:
                        f.write(f"Tracking target(Galatic): Az = {self.point_az}; Alt = {self.point_alt}. l = {input_l}; b = {input_b}.\n ")
                        self.file_write_count += 1

            else:
                print("......")
        else:
            if self.ui.obj_name.currentText() == "Sun":
                self.rc.goto_sun()
                if self.file_write_count == 0:
                    with open(file_name,"a") as f:
                        f.write(f"Tracking target: Sun.\n")
                        self.file_write_count += 1

            elif self.ui.obj_name.currentText() == "Moon":
                self.rc.goto_moon()
                if self.file_write_count == 0:
                    with open(file_name,"a") as f:
                        f.write(f"Tracking target: Sun.\n")
                        self.file_write_count += 1

            elif self.ui.obj_name.currentText() == "M31":
                print("Goto M31.....")
            else:
                print("GoTo by coords...")

    def tracking_timer(self):
        self.file_write_count = 0
        self.lock_radec = False
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f"Tracking start time:{temp_time}.\n")
            f.write("**********\n")
        self.tracking()
        t_s = int( float( self.ui.track_time.toPlainText() )*3600*1000 )
        self.tracking_time_len = QtCore.QTimer(self)
        self.tracking_time_len.singleShot(t_s,self.tracking_timer_end)


    def tracking(self):
        print("Rotator start tracking target.....\n")
        self.timer_track = QtCore.QTimer(self)
        self.timer_track.timeout.connect(self.tracking_obj)
        self.timer_track.start(1000)

    def tracking_timer_end(self):
        temp_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f"Tracking stop time:{temp_time}.\n")
            f.write("**********\n")
        self.timer_track.stop()
        print("Rotator stopped.\n")
        

    # start timer
    def timer_start(self):
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.plot_sky)
        self.timer.start(1000)

    def timer_end(self):
        self.timer.stop()
        ctime1 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        ctime = Time(ctime1)  # Eastern Daylight Time 
        print(f"Stop at {ctime}.\n\n")

        with open(file_name,"a") as f:
            f.write("**********\n")
            f.write(f"Tracking start time:{ctime1}. ")
            f.write("**********\n")

    def plot_sky(self):

        ctime1,xl,yl,lr,br,sun_az,sun_alt,moon_az,moon_alt = set_calculate_sky(self.loc_alt,self.loc_lon,self.loc_hei,self.loc_name,self.loc_tiz)

        self.ui.sky_area.clear()
        self.ui.sky_area.setLimits(xMin=-190,xMax=190,yMin=0,yMax=95)
        self.ui.sky_area.setXRange(min=-190,max=190)
        self.ui.sky_area.setYRange(min=0,max=95)

        sun = pg.TextItem("Sun",color='#f28500') #"#000000")
        self.ui.sky_area.addItem(sun)
        sun.setPos(sun_az+3, sun_alt+3)

        moon = pg.TextItem("Moon",color="#778899")
        self.ui.sky_area.addItem(moon)
        moon.setPos(moon_az+3, moon_alt+3)

        for i in range(len(yl)):
            name = pg.TextItem(f"({ int(br[i]) },{ lr[i] })",color="#00bfff")
            self.ui.sky_area.addItem(name)
            name.setPos(xl[i]+3, yl[i]+3)

        self.ui.sky_area.plot().setData([sun_az],[sun_alt],pen='r', symbol='o',symbolSize=15,symbolBrush=(218,165,32))
        self.ui.sky_area.plot().setData([moon_az],[moon_alt],pen='r', symbol='o',symbolSize=15,symbolBrush=(245,245,245))
        self.ui.sky_area.plot().setData(xl,yl,pen=None, symbol='s',symbolSize=8,symbolBrush=(65,105,255))
        self.ui.loc_info.setText("Location: "+self.loc_name+"\n\n"+"Current Time: "+str(ctime1))

        if self.ui.obj_name.currentText() == "Sun":

            self.ui.sky_area.plot().setData([sun_az],[sun_alt],pen='r', symbol='+',symbolSize=15,symbolBrush=(128,0,0))

            telescope = pg.TextItem("Telescope",color="#800000")
            self.ui.sky_area.addItem(telescope)
            telescope.setPos(sun_az+2, sun_alt+5)


        elif self.ui.obj_name.currentText() == "Moon":
            self.ui.sky_area.plot().setData([moon_az],[moon_alt],pen='r', symbol='+',symbolSize=15,symbolBrush=(128,0,0))

            telescope = pg.TextItem("Telescope",color="#800000")
            self.ui.sky_area.addItem(telescope)
            telescope.setPos(moon_az+2, moon_alt+5)

        else:
            self.ui.sky_area.plot().setData([self.point_az-180],[self.point_alt],pen='r', symbol='+',symbolSize=15,symbolBrush=(128,0,0))
            telescope = pg.TextItem("Telescope",color="#800000")
            self.ui.sky_area.addItem(telescope)
            telescope.setPos(self.point_az-180+2, self.point_alt+5)
        #cal_sky.join()

if __name__ == '__main__':
    app_set = QtWidgets.QApplication([sys.argv])
    mainw_set = setting_location()
    mainw_set.show()
    app_set.exec_()
    time.sleep(0.2)

    app = QtWidgets.QApplication([sys.argv])
    mainw = Rc_Window(mainw_set.loc,mainw_set.alt,mainw_set.lon,mainw_set.hei,mainw_set.tiz,mainw_set.port)
    mainw.show()
    sys.exit(app.exec_())
