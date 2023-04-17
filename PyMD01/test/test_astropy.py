#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:54:01 2023

@author: jiomer
"""

#%%
import os
import sys
import time
import serial
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


#%%

Boston = EarthLocation(lat=42.3*u.deg, lon=-71*u.deg, height=0*u.m) 

utcoffset = -5*u.hour
ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset

#%%


obj_rd = SkyCoord('03h46m30s', '+23d56m57s', frame='icrs')
print(f"Your input coordinates is: \n  RA:{obj_rd.ra}. \n  DEC:{obj_rd.dec}.")


altaz = AltAz(obstime=ctime, location=Boston)
obj_ae = obj_rd.transform_to(altaz)
print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")

#%%

obj_rd = SkyCoord('03h46m30s', '+23d56m57s', frame='icrs')
print(f"Your input coordinates is: \n  RA:{obj_rd.ra}. \n  DEC:{obj_rd.dec}.")
utcoffset = -5*u.hour
ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset

altaz = AltAz(obstime=ctime, location=Boston)
obj_ae = obj_rd.transform_to(altaz)
print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")


#%%

from threading import Thread
import time


#%%
def jisuan(a,b):
    i = 0
    count = 5
    while i <= count:
        a = b/a
        b = a/b
        print(f"This is {i} times calculation, resultis a = {a}, b = {b}.\n")
        i += 1
        time.sleep(2)
        
        
def print_num(n1,n2):
    j = 0
    while j <= 6:
        time.sleep(1)
        print(f"Here is another thread, the number are: {n1} and {n2}...\n")
        j += 1
        
#%%
#if __name__ == "__main__":
num1 = 2
num2 = 3
j = 0
print("This is main thread...\nWaiting the threadings...\n")
t1 = Thread(target=jisuan, args=(num1,num2))
t2 = Thread(target=print_num, args=(num1,num2))
    
t1.start()
t2.start()
    
t1.join()
t2.join()

print("This is still main thread...")







#%%
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun, get_moon
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

#%%
bos = EarthLocation(lon=-71.06 * u.deg, lat=42.36 * u.deg, height=1400 * u.m)
midnight = Time('2023-01-09 04:00:00') - (-5*u.hour)


delta_midnight = np.linspace(0, 24, 240)*u.hour
times_July12_to_13 = midnight + delta_midnight
frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bos)
sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)

fig = plt.figure(figsize=(14,12))
plt.scatter( sunaltazs_July12_to_13.az.value, sunaltazs_July12_to_13.alt.value )


#%%

bos = EarthLocation(lon=-71.06 * u.deg, lat=42.36 * u.deg, height=1400 * u.m)
midnight = Time('2023-01-09 16:17:00') - (-5*u.hour)


delta_midnight = np.linspace(0, 90, 1200)*u.hour
dates = midnight + delta_midnight
moon_dates = AltAz(obstime=dates, location=bos)
moon_loc = get_moon(dates).transform_to(moon_dates)

fig = plt.figure(figsize=(14,12))
plt.scatter( moon_loc.az.value, moon_loc.alt.value )


#%%
import sys

from PySide2.QtWidgets import QApplication, QDialog, QMainWindow, QPushButton


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("My App")

        button = QPushButton("Press me for a dialog!")
        button.clicked.connect(self.button_clicked)
        self.setCentralWidget(button)

    def button_clicked(self, s):
        print("click", s)

        dlg = QDialog(self)
        dlg.setWindowTitle("HELLO!")
        dlg.exec_()


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec_()









































































# %%
