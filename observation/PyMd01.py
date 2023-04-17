#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 13:38:58 2023

@author: jiomer

Using Python control the SPID MD-01 rotator.

"""


import os
import sys
import time
import serial
import argparse
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic
from astropy.coordinates import get_sun, get_moon
import numpy as np
#from threading import Thread


lon = -71.06
lat = 43.36
height = 30
timezone = -4
portname = '/dev/ttyUSB0'



def num_range(num,lower,upper):
    if num < lower:
        r_num = lower
        #print(f"The number you input {num} is less than lower range, re-set number equals to {lower}.")
    elif num > upper:
        r_num = upper
        #print(f"The number you input {num} is larger than upper range, re-set number equals to {upper}.")
    else:
        r_num = num
        #print("Now set the number is {num}.")
        
    return r_num





five_0 = chr(0) + chr(0) + chr(0) + chr(0) + chr(0)
four_0 = chr(0) + chr(0) + chr(0) + chr(0)

require_up = chr(0x57) + chr(0x04) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_down = chr(0x57) + chr(0x08) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_left = chr(0x57) + chr(0x01) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_right = chr(0x57) + chr(0x02) + five_0 + four_0 + chr(0x14) + chr(0x20)

require_leftup = chr(0x57) + chr(0x05) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_leftdown = chr(0x57) + chr(0x09) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_rightup = chr(0x57) + chr(0x06) + five_0 + four_0 + chr(0x14) + chr(0x20)
require_rightdown = chr(0x57) + chr(0x0A) + five_0 + four_0 + chr(0x14) + chr(0x20)


class RotatorController:
    """
    This class is working for commu with the rotator and control it pointing.
    And then track the sky objet.
    
    """
    def __init__(self,port: str = portname, baud=115200, wait_time=1,loclon=lon,loclat=lat,locheight=height,utcoffset=timezone,daz=0,dalt=0):
        """
        port: Serial connected
        baud: Rotator manual suggested equals to 115200
        wait_time: set a timeout value
        response: rotator response
        az: rotator current azimuth
        el: rotator current elevation
        daz: delta angle in az, from north to east is +
        dalt: delta angle in alt, from 0 to zenith is +
        
        """
        
        self.port = port
        try:
            os.stat(self.port)
        except OSError as oserror:
            raise ValueError(f"Given {self.port}, which does not exist!") from oserror
            sys.exit(1)
            
        self.baud = baud
        self.wait_time = wait_time
        self.ser = serial.Serial(self.port, self.baud, timeout = self.wait_time)

        self.response = self.get_response()
        self.az_multi = self.response[5]
        self.el_multi = self.response[10]
        #self.response_100 = self.get_response_100()
        
        self.angle_res_10 = 0.1 # degree
        self.angle_res_100 = 0.01 # degree
        
        self.location = EarthLocation(lon=loclon * u.deg, lat=loclat * u.deg, height=locheight * u.m)
        self.utcoffset = utcoffset

        self.daz = daz
        self.dalt = dalt
        
    def get_response(self,try_time=10):
        """
        Try to get response from the rotator.
        You should do it 10 times.
        Position resolution is 0.01
        
        """
        #print(f"Stat of the serial is {self.ser.isOpen()}.")
        try_time = int(try_time)
        
        require = chr(0x57) + five_0 + five_0 + chr(0x1f) + chr(0x20)
        
        count = 0
        response = b""
        while count < try_time:
            count += 1
            out_stat = self.ser.write(require.encode())
            time.sleep(0.5)
            response = self.ser.read(12)
            if len(response) >= 12:
                #print(out_stat)
                break
        else:
            print(f"Try {try_time} times but no response, ended.")
        return response
    
    def get_response_100(self,try_time=10):
        """
        Try to get response from the rotator.
        You should do it 10 times.
        Position resolution is 0.01
        
        """
        #print(f"Stat of the serial is {self.ser.isOpen()}.")
        try_time = int(try_time)
        
        require = chr(0x57) + five_0 + five_0 + chr(0x6f) + chr(0x20)
        
        count = 0
        response_100 = b""
        while count < try_time:
            count += 1
            out_stat = self.ser.write(require.encode())
            time.sleep(0.5)
            response_100 = self.ser.read(12)
            if len(response_100) >= 12:
                #print(out_stat)
                break
        else:
            print(f"Try {try_time} times but no response, ended.")
        return response_100
    
    def get_loc(self):
        """
        From response get rotator position
        
        """
        time.sleep(0.1)
        i = 0
        while i < 5:
            self.response = self.get_response()
            self.az = self.response[1]*100+self.response[2]*10 + self.response[3] + self.response[4]/10 - 360.0 + self.daz
            self.el = self.response[6]*100+self.response[7]*10 + self.response[8] + self.response[9]/10 - 360.0 + self.dalt
            time.sleep(0.1/5)
            i += 1
        print(f"rotator current position is AZ = {self.az}, EL = {self.el}.")
        return self.az - self.daz, self.el - self.dalt
    
    
    def get_loc_100(self):
        """
        From response get rotator position, resolution = 0.01 degree.
        
        """
        time.sleep(1)
        i = 0
        while i < 5:
            self.response_100 = self.get_response_100()
            self.az_100 = self.response_100[1]*100+self.response_100[2]*10 + self.response_100[3] + self.response_100[4]/10 + self.response_100[5]/100 - 360.0 + self.daz
            self.el_100 = self.response_100[6]*100+self.response_100[7]*10 + self.response_100[8] + self.response_100[9]/10 + self.response_100[10]/100 - 360.0 + self.dalt
            time.sleep(0.1/5)
            i += 1
        print(f"rotator current position is AZ = {self.az_100}, EL = {self.el_100}.")
        return self.az_100 - self.daz, self.el_100 - self.dalt
    
    
    def set_power(self,az_power:int=100,el_power:int=100):
        """
        Set motor power.
        
        """
        az_power = num_range(az_power, 0, 100)
        el_power = num_range(el_power, 0, 100)
        four_0 = chr(0) + chr(0) + chr(0) + chr(0)
        power_set = chr(0x57) + four_0 + chr( int(az_power) ) + four_0 + chr( int(el_power) ) + chr(0xf7) + chr(0x20)
        out_stat = self.ser.write(power_set.encode())
        #print(out_stat)
        time.sleep(0.2)
        print(f"Now the az motor power is {int(az_power)}%, and el motor power is {int(el_power)}%.")
        return True
    
    
    
    def AZEL_goto(self, az: float = None, el: float = None) -> bool:
        """
        Communicates with the rotator to set the desired pointing location.

        args:
            az - desired azimuth, must be between [0, 180]

            el - desired elivation
        returns
            success - weather the message was sent.
        """
        if not 0 <= el <= 180:
            return False
        if not 0 <= az <= 360:
            return False

        az = az + 360 - self.daz
        el = 180 - el - self.dalt + 360

        az *= self.az_multi
        el *= self.el_multi

        az = str(int(az))
        el = str(int(el))

        if len(az) == 3:
            az = "0" + az
        if len(el) == 3:
            el = "0" + el
        #print(az,el)
        # Message to send
        message = (
            chr(0x57)
            + az
            + chr(self.az_multi)
            + el
            + chr(self.el_multi)
            + chr(47)
            + chr(32)
        )

        # Send message
        out_stat = self.ser.write(message.encode())
        return True

    
    def RaDec_goto(self, RA:str=None, DEC:str=None):
        """
        Go to a coordinate in  RA, Dec
        obj_rd: sky object coordinates ra dec.
        obj_ae: sky object coordinates az el.
        az in range (0,360)
        el in range (0,180)
        
        """
        obj_rd = SkyCoord(RA, DEC, frame='icrs')
        #print(f"Your input coordinates is: \n  RA:{obj_rd.ra}. \n  DEC:{obj_rd.dec}.")
        # time area, west is - and east is +
        # for example Boston at UTC -5h so that utcoffset = -5h
        utcoffset = self.utcoffset*u.hour
        ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset
        
        altaz = AltAz(obstime=ctime, location=self.location)
        obj_ae = obj_rd.transform_to(altaz)
        #print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")
        
        if 0<= obj_ae.alt.degree < 180:
            self.AZEL_goto(obj_ae.az.degree,obj_ae.alt.degree)
        
        else:
            print("The sky object is below the horizon ..... \nPlz change to another object.")
        
        return obj_ae.az.value,obj_ae.alt.value
    
    def galatic_goto(self, L:str=None, B:str=None):
        # 将银道坐标转换为天球坐标
        c = SkyCoord(l=L, b=B, frame=Galactic, unit=(u.deg, u.deg)).icrs
        RA_goto = c.ra.to_string(unit=u.hourangle, sep='hms', precision=1, pad=True)
        DEC_goto = c.dec.to_string(unit=u.degree, sep='dms', precision=1, pad=True)
        self.RaDec_goto(RA=RA_goto,DEC=DEC_goto)
        return True
    


    def AZEL_goto_100(self, az: float = None, el: float = None) -> bool:
        """
        Communicates with the rotator to set the desired pointing location.

        args:
            az - desired azimuth, must be between [0, 180]

            el - desired elivation
        returns
            success - weather the message was sent.
        """
        az = az + 360 - self.daz
        el = 180 - el - self.dalt + 360

        az *= 100
        el *= 100

        az = str(int(az))
        el = str(int(el))

        message = (
            chr(0x57)
            + az
            + el
            + chr(0x5f)
            + chr(0x20)
        )

        # Send message
        out_stat = self.ser.write(message.encode())
        return True
    
    
    def goto_moon(self):
        """
        Tracking moon.
        
        """
        
        utcoffset = self.utcoffset*u.hour
        ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset
        loc_moon_azel = AltAz(obstime=ctime, location=self.location)
        moon_azel = get_moon(ctime).transform_to(loc_moon_azel)
        
        if 0<= moon_azel.alt.value < 180:
            self.AZEL_goto_100(moon_azel.az.value,moon_azel.alt.value)
        
        else:
            self.stop()
            print("The Sun is below the horizon now ......")
        
        return True

    def tracking_moon(self,hours:float = 0.0):
        """
        Tracking the sun position.
        hours = 0 is tracking without time limits.
        
        """
        if hours == 0.0:
            utcoffset = self.utcoffset*u.hour
            while True:
                ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset
                loc_moon_azel = AltAz(obstime=ctime, location=self.location)
                moon_azel = get_moon(ctime).transform_to(loc_moon_azel)
                
                if 0<= moon_azel.alt.value < 180:
                    self.AZEL_goto_100(moon_azel.az.value,moon_azel.alt.value)
                    time.sleep(5)
                else:
                    break
                    self.stop()
                    print("The Moon is below the horizon now ......")
                    
        else:
            ftime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) + (hours - self.utcoffset)*u.hour
            ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
            while ctime < ftime:
                ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
                loc_sun_azel = AltAz(obstime=ctime, location=self.location)
                sun_azel = get_moon(ctime).transform_to(loc_sun_azel)
                
                if 0<= sun_azel.alt.value < 180:
                    self.AZEL_goto_100(sun_azel.az.value,sun_azel.alt.value)
                    time.sleep(5)
                else:
                    break
                    self.stop()
                    print("The Sun is below the horizon now ......")
                    
            self.stop()
                
        return True



    def goto_sun(self):
        """
        goto position of sun.

        """
        utcoffset = self.utcoffset*u.hour
        ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset
        loc_sun_azel = AltAz(obstime=ctime, location=self.location)
        sun_azel = get_sun(ctime).transform_to(loc_sun_azel)
        
        if 0<= sun_azel.alt.value < 180:
            self.AZEL_goto_100(sun_azel.az.value,sun_azel.alt.value)
        
        else:
            self.stop()
            print("The Sun is below the horizon now ......")
        
        return True
    
    
    def tracking_sun(self,hours: float = 0.0):
        """
        Tracking the sun position.
        hours = 0 is tracking without time limits.
        
        """
        if hours == 0.0:
            utcoffset = self.utcoffset*u.hour
            while True:
                ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset
                loc_sun_azel = AltAz(obstime=ctime, location=self.location)
                sun_azel = get_sun(ctime).transform_to(loc_sun_azel)
                
                if 0<= sun_azel.alt.value < 180:
                    self.AZEL_goto_100(sun_azel.az.value,sun_azel.alt.value)
                    time.sleep(5)
                else:
                    break
                    self.stop()
                    print("The Sun is below the horizon now ......")
                    
        else:
            ftime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) + (hours - self.utcoffset)*u.hour
            ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
            while ctime < ftime:
                ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
                loc_sun_azel = AltAz(obstime=ctime, location=self.location)
                sun_azel = get_sun(ctime).transform_to(loc_sun_azel)
                
                if 0<= sun_azel.alt.value < 180:
                    self.AZEL_goto_100(sun_azel.az.value,sun_azel.alt.value)
                    time.sleep(5)
                else:
                    break
                    self.stop()
                    print("The Sun is below the horizon now ......")
                    
            self.stop()
                
        return True
    
    
    def tracking(self,RA:str=None,DEC:str=None,hours:float=None):
        """
        虽然rotator可以以0.01度的精度运转，但是连续使用0.01度的分辨率驱动rotator时每经过0.1度才会运动。
        不知道能不能改。
        
        """
        print("Now start tracking......\n")
        f_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) + (hours - self.utcoffset)*u.hour
        current_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
        obj_rd = SkyCoord(RA, DEC, frame='icrs')
        while f_time > current_time:
            current_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
            altaz = AltAz(obstime=current_time, location=self.location)
            obj_ae = obj_rd.transform_to(altaz)
            #print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")
                
            if 0<= obj_ae.alt.degree < 180:
                self.AZEL_goto_100(obj_ae.az.degree,obj_ae.alt.degree)
                time.sleep(5)
                
            else:
                self.stop()
                print("The sky object is below the horizon ..... \n")
                break

        self.stop()
        print(f"Total tracking period is {hours} hours. Now finished.")
        
        return True
    
    def tracking_galatic(self,L:str=None,B:str=None,hours:float=None):
        print("Now start tracking in galatic coordinates......\n")
        f_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) + (hours - self.utcoffset)*u.hour
        current_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
        # 将银道坐标转换为天球坐标
        c = SkyCoord(l=L, b=B, frame=Galactic, unit=(u.deg, u.deg)).icrs
        RA = c.ra.to_string(unit=u.hourangle, sep='hms', precision=1, pad=True)
        DEC = c.dec.to_string(unit=u.degree, sep='dms', precision=1, pad=True)

        obj_rd = SkyCoord(RA, DEC, frame='icrs')
        while f_time > current_time:
            current_time = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - self.utcoffset*u.hour
            altaz = AltAz(obstime=current_time, location=self.location)
            obj_ae = obj_rd.transform_to(altaz)
            #print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")
                
            if 0<= obj_ae.alt.degree < 180:
                self.AZEL_goto_100(obj_ae.az.degree,obj_ae.alt.degree)
                time.sleep(5)
                
            else:
                self.stop()
                print("The sky object is below the horizon ..... \n")
                break

        self.stop()
        print(f"Total tracking period is {hours} hours. Now finished.")
        return True
    
    def longtime_simulation(self, RA:str=None, DEC:str=None,obslength:float=6,sleep_time:float=0.5):
        """
        Go to a coordinate in  RA, Dec
        obj_rd: sky object coordinates ra dec.
        obj_ae: sky object coordinates az el.
        az in range (0,360)
        el in range (0,180)
        
        """
        obj_rd = SkyCoord(RA, DEC, frame='icrs')
        #print(f"Your input coordinates is: \n  RA:{obj_rd.ra}. \n  DEC:{obj_rd.dec}.")
        # time area, west is - and east is +
        # for example Boston at UTC -5h so that utcoffset = -5h
        utcoffset = self.utcoffset*u.hour
        ctime = Time( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) - utcoffset + np.arange(0,obslength,0.01)*u.h
        
        altaz = AltAz(obstime=ctime, location=self.location)
        obj_ae = obj_rd.transform_to(altaz)
        #print(f"Now the AZEL coordinates is: \n  AZ:{obj_ae.az}. \n  EL:{obj_ae.alt}.")
        
        for i in range(len(ctime)):
            if 0<= obj_ae.alt.degree < 180:
                self.AZEL_goto(obj_ae[i].az.degree,obj_ae[i].alt.degree)
                time.sleep(sleep_time)
        
            else:
                print("The sky object is below the horizon ..... \nPlz change to another object.")
                break


    def move_up(self):
        print("move up")
        out_stat = self.ser.write(require_up.encode())

    def move_left_up(self):
        print("move left up")
        out_stat = self.ser.write(require_leftup.encode())

    def move_left(self):
        print("move left")
        out_stat = self.ser.write(require_left.encode())

    def move_left_down(self):
        print("move left down")
        out_stat = self.ser.write(require_leftdown.encode())

    def move_down(self):
        print("move down")
        out_stat = self.ser.write(require_down.encode())

    def move_right_down(self):
        print("move right down")
        out_stat = self.ser.write(require_rightdown.encode())

    def move_right(self):
        print("move right")
        out_stat = self.ser.write(require_right.encode())

    def move_right_up(self):
        print("move right up")
        out_stat = self.ser.write(require_rightup.encode())


    def stop(self):
        """
        Stops the rotator.
        
        """
        stop_message = chr(87) + five_0 + five_0 + chr(15) + chr(32)
        out_stat = self.ser.write(stop_message.encode())
        time.sleep(0.1)
        print("Rotator stopped.")
        
    
    def close_comm(self,return_home=False):
        if return_home == True:
            self.AZEL_goto(0,0)
        else:
            pass
        time.sleep(10)
        self.ser.close()
        print(f"Now {self.port} is closed.")
        #return True
























