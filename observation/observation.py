
#%%
import os
import time
import subprocess
import numpy as np
import threading
from datetime import datetime, timedelta

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic
import astropy.units as u
from astropy.time import Time

import configparser

from PyMd01 import RotatorController

#%%

def run_sdr(sdr_freq,sdr_bw,sdr_g,sdr_o,sdr_fs,sdr_fn1,sdr_fn2,sdr_path):
    subprocess.run([
        '/usr/bin/python3','./RTLSDR2.py','-f',str(sdr_freq),
        '-b',str(sdr_bw),'-g',str(sdr_g),'-o',str(sdr_o), '-fs',
        str(sdr_fs), '-n1',sdr_fn1,'-n2',sdr_fn2, '-p',sdr_path
                    ])


def radec_dist(ra1,dec1,ra2,dec2,inumber):
    bos = EarthLocation(lat=42.3601*u.deg, lon=-71.0942*u.deg, height=30*u.m)
    bos_time = datetime.now() - timedelta(hours=-4)
    coord1 = SkyCoord(ra1,dec1,frame='icrs')
    coord2 = SkyCoord(ra2,dec2,frame='icrs')
    
    # 将天赤道坐标转换为地平坐标系
    altaz1 = coord1.transform_to(AltAz(obstime=bos_time, location=bos))
    altaz2 = coord2.transform_to(AltAz(obstime=bos_time, location=bos))
    #print(altaz1,altaz2)
    delta_az = max(altaz1.az.deg, altaz2.az.deg) - min(altaz1.az.deg, altaz2.az.deg)
    delta_alt = max(altaz1.alt.deg, altaz2.alt.deg) - min(altaz1.alt.deg, altaz2.alt.deg)
    #print(delta_az,delta_alt)

    max_angle_dist = np.max( [delta_az,delta_alt] )
    #print(max_angle_dist)
    spend_time = int( np.round( max_angle_dist/1.5 ) )
    if inumber == 0:
        spend_time += 22
    else:
        pass

    if spend_time == 0:
        spend_time = 1
    else:
        pass

    if spend_time > int(360/1.5)+1:
        spend_time = int(360/1.5)+1
    else:
        pass

    return spend_time


def grid_coords(c_ra,c_dec,n_point,gap_deg):
    """
    c_ra: center point ra
    c_dec: center point dec
    n_point: points number for a side, total points is n_point*n_point
    gap_deg: distance in degree between two point
    onetime: observation time in one point.
    """
    grid_ra_list = []
    grid_dec_list = []

    center_coord= SkyCoord(ra=c_ra,dec=c_dec)
    c_ra_deg = center_coord.ra.value
    c_dec_deg = center_coord.dec.value
    temp_ralist = np.arange( -(n_point-1)/2,(n_point-1)/2+1,1 )*gap_deg
    temp_declist = []
    for i in range(n_point):
        temp_decl = []
        temp_decl += [i+-(n_point-1)/2]*n_point
        temp_declist.extend(temp_decl)
    temp_declist = np.array(temp_declist)*gap_deg

    ra_list = c_ra_deg + temp_ralist
    ra_list = np.tile(ra_list,n_point)
    dec_list = temp_declist + c_dec_deg

    grid_coords = SkyCoord(ra=ra_list*u.deg,dec=dec_list*u.deg)
    grid_str_coords = grid_coords.to_string('hmsdms')

    for i in range(len(grid_str_coords)):
        grid_ra_list.append(grid_str_coords[i][:9])
        grid_dec_list.append(grid_str_coords[i][10:])

    #print(grid_ra_list,grid_dec_list)

    return grid_ra_list,grid_dec_list

#########################################################################
#########################################################################
#########################################################################

#%%

# single observation

temp_date = datetime.now()
Date = temp_date.strftime('%Y-%m-%d')
Time = temp_date.strftime('%H:%M:%S')
cfreq = 1420e6
bandwidth = 1.8e6
gain = 50
obstime = 300
fftsize = 4096
path = './data/single/'
fn1bl = 'SDR1-bl.txt'
fn2bl = 'SDR2-bl.txt'
fn1 = 'SDR1-data.txt'
fn2 = 'SDR2-data.txt'
coord_format = 'EQ: ra_dec'

log_name = Date + "_SDR_single_point.ini"

start_ra = '0h0m0s'
start_dec = '90d0m0s'

single_blra = '8h0m12s'
single_bldec = '73d55m20s'

single_ra = '7h34m35s'
single_dec = '-27d01m18s'

#%%

rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)

#%%
# 创建配置文件对象
config = configparser.ConfigParser()

# sdr_freq,sdr_bw,sdr_g,sdr_o,sdr_fs,sdr_fn1,sdr_fn2,sdr_path
t_bl = threading.Thread( target=run_sdr,
args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1bl,fn2bl,path,) )


#################
# record baseline
#################
rc.RaDec_goto(single_blra,single_bldec)
sleeptime = radec_dist(start_ra,start_dec,single_blra,single_bldec,1)
time.sleep(sleeptime)

temp_date = datetime.now()
obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
t_bl.start()
rc.tracking(single_blra,single_bldec,obstime/3600)
t_bl.join()


# 添加一个节和选项
config['SDR1-bl_'+ Time ] = {
'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn1bl,
'observation_mode': 'single_point', 'observation_time': obser_time, 
'coord_format': coord_format, 'coords1': single_blra, 'coords2': single_bldec
}

config['SDR2-bl_' + Time ] = {
'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn2bl,
'observation_mode': 'single_point', 'observation_time': obser_time,
'coord_format': coord_format, 'coords1': single_blra, 'coords2': single_bldec
}

# 将配置写入到文件中
with open(log_name, 'a') as configfile:
    config.write(configfile)

time.sleep(0.1)

#%%

#################
# data record
#################
t_data = threading.Thread( target=run_sdr,
args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )

config = configparser.ConfigParser()
# record data
rc.RaDec_goto(single_ra,single_dec)
sleeptime = radec_dist(single_blra,single_bldec,single_ra,single_dec,1)
time.sleep(sleeptime)

temp_date = datetime.now()
obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
t_data.start()
rc.tracking(single_ra,single_dec,obstime/3600)
t_data.join()

# 添加一个节和选项
config['SDR1-data_' + Time ] = {
'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
'observation_mode': 'single_point', 'observation_time': obser_time,
'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
}

config['SDR2-data_' + Time ] = {
'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
'observation_mode': 'single_point', 'observation_time': obser_time,
'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
}

# 将配置写入到文件中
with open(log_name, 'a') as configfile:
    config.write(configfile)

#%%

rc.close_comm()

# single observation end.

############################################################################
############################################################################
############################################################################

#%%

# survey whole sky

start_ra = '0h0m0s'
start_dec = '90d0m0s'

temp_date = datetime.now()
Date = temp_date.strftime('%Y-%m-%d')
#Time = temp_date.strftime('%H:%M:%S')
cfreq = 1420e6
bandwidth = 1.8e6
gain = 50
obstime = 30
fftsize = 4096
path = './data/survey/'
fn1bl = 'SDR1-survey_bl.txt'
fn2bl = 'SDR2-survey_bl.txt'
fn1 = 'SDR1-survey_data.txt'
fn2 = 'SDR2-survey_data.txt'
coord_format = 'EQ: ra_dec'

log_name = Date + "_SDR_survey.ini"

# number of how many data collected to goto record baseline
bl_gap = 5

blra = '6h0m0s'
bldec = '26d0m0s'

#%%

subprocess.run(['python', './survey_plan.py'])

#%%

# 获取当前文件夹路径
folder_path = os.getcwd()

# 获取当前文件夹下所有以.npz为扩展名的文件
files = [file for file in os.listdir(folder_path) if file.endswith('.npz')]

# 读取每个文件并将其存储在字典中
data_dict = {}
for i in range(len(files)):
    file_path = os.path.join(folder_path, files[i])
    data_dict[i] = np.load(file_path)

#print(data_dict)

obs_ra = data_dict[0]["ra"]
obs_dec = data_dict[0]["dec"]
obs_settime = data_dict[0]["settime"]

#%%

obs_str_ra = []
obs_str_dec = []

# create SkyCoord object
coords = SkyCoord(ra=obs_ra*u.deg, dec=obs_dec*u.deg, frame='icrs')

# 转换为hmsdms格式的字符串，并取出ra和dec
for c in coords:
    temp_str = c.to_string('hmsdms')
    ra, dec = temp_str.split()[0], temp_str.split()[1]
    obs_str_ra.append(ra)
    obs_str_dec.append(dec)

obs_settime = np.array(obs_settime)
obs_str_ra = np.array(obs_str_ra)
obs_str_dec = np.array(obs_str_dec)

# 获取array3数组的升序排列的索引值
sorted_indices = np.argsort(obs_settime)

# 使用sorted_indices对array1和array2进行重新排序
obs_settime = obs_settime[sorted_indices]
obs_str_ra = obs_str_ra[sorted_indices]
obs_str_dec = obs_str_dec[sorted_indices]


#%%

rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)


#%%

for i in range( 5 ):
    if i % bl_gap == 0:
        print(i)
        config = configparser.ConfigParser()
        # sdr_freq,sdr_bw,sdr_g,sdr_o,sdr_fs,sdr_fn1,sdr_fn2,sdr_path
        t_bl = threading.Thread( target=run_sdr,
        args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),
              fn1bl,fn2bl,path,) )
        
        sleeptime = radec_dist( start_ra,start_dec,blra,bldec,i )
        print(sleeptime)
        
        # record baseline
        rc.RaDec_goto(blra,bldec)
        time.sleep(sleeptime)

        temp_date = datetime.now()
        obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
        t_bl.start()
        rc.tracking(blra,bldec,obstime/3600)
        t_bl.join()

        start_ra = blra
        start_dec = bldec


        # 添加一个节和选项
        config['SDR1-bl_'+ str(i) ] = {
        'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
        'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn1bl,
        'observation_mode': 'survey', 'observation_time': obser_time, 
        'coord_format': coord_format, 'coords1': blra, 'coords2': bldec
        }

        config['SDR2-bl_' + str(i) ] = {
        'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
        'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn2bl,
        'observation_mode': 'survey', 'observation_time': obser_time,
        'coord_format': coord_format, 'coords1': blra, 'coords2': bldec
        }

        # 将配置写入到文件中
        with open(log_name, 'a') as configfile:
            config.write(configfile)
    else:
        pass
    
    single_ra = obs_str_ra[i]
    single_dec = obs_str_dec[i]

    sleeptime = radec_dist(start_ra,start_dec,single_ra,single_dec,i)
    print(sleeptime)

    start_ra = single_ra
    start_dec = single_dec

    t_data = threading.Thread( target=run_sdr,
    args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )

    config = configparser.ConfigParser()
    # record data
    rc.RaDec_goto(single_ra,single_dec)
    time.sleep(sleeptime)

    temp_date = datetime.now()
    obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
    t_data.start()
    rc.tracking(single_ra,single_dec,obstime/3600)
    t_data.join()

    # 添加一个节和选项
    config['SDR1-data_' + str(i) ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
    'observation_mode': 'survey', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
    }

    config['SDR2-data_' + str(i) ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
    'observation_mode': 'survey', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
    }

    # 将配置写入到文件中
    with open(log_name, 'a') as configfile:
        config.write(configfile)

# %%
rc.close_comm()


###########################################################################
###########################################################################
###########################################################################


# %%

# grid observation

temp_date = datetime.now()
Date = temp_date.strftime('%Y-%m-%d')
Time = temp_date.strftime('%H:%M:%S')
cfreq = 1420e6
bandwidth = 1.8e6
npoints=5
gapdeg=7
gain = 50
obstime = 30
fftsize = 4096
path = './data/grid/'

fn1 = 'SDR1-grid_data.txt'
fn2 = 'SDR2-grid_data.txt'
coord_format = 'EQ: ra_dec'

log_name = Date + "_SDR_grid_points.ini"

start_ra = '1h00m54s'
start_dec = '7d12m43s'

center_ra = '1h7m55s'
center_dec = '7d12m52s'

#%%

rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)

#%%

# record data
rc.RaDec_goto(center_ra,center_dec)
sleeptime = radec_dist(center_ra,center_dec,start_ra,start_dec,0)
print(sleeptime)
#time.sleep(8)
time.sleep(sleeptime)

grids_ra,grids_dec = grid_coords(center_ra,center_dec,npoints,gapdeg)


for i in range(len(grids_ra)):
    config = configparser.ConfigParser()
    temp_date = datetime.now()
    obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')

    rc.RaDec_goto( grids_ra[i],grids_dec[i] )
    time.sleep( int( np.round(gapdeg/1.5)+1 ) )
    t_data = threading.Thread( target=run_sdr,
    args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )
    t_data.start()
    rc.tracking(grids_ra[i],grids_dec[i],obstime/3600)
    t_data.join()

    # 添加一个节和选项
    config['SDR1-data_'+ str(i+1) + "_" + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': grids_ra[i], 'coords2': grids_dec[i]
    }

    config['SDR2-data_' + str(i+1) + "_" + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': grids_ra[i], 'coords2': grids_dec[i]
    }

    # 将配置写入到文件中
    with open(log_name, 'a') as configfile:
        config.write(configfile)



# %%
def obs_bl(cfreq,bandwidth,gain,obstime,fftsize,path,fn1bl,fn2bl,start_ra,start_dec,single_blra,single_bldec):
    temp_date = datetime.now()
    Date = temp_date.strftime('%Y-%m-%d')
    Time = temp_date.strftime('%H:%M:%S')
    #cfreq = 1420.4e6
    #bandwidth = 1.8e6
    #gain = 50
    #obstime = 300
    #fftsize = 4096
    #path = './data/single/'
    #fn1bl = 'SDR1-bl.txt'
    #fn2bl = 'SDR2-bl.txt'
    coord_format = 'EQ: ra_dec'
    log_name = Date + "_SDR_single_point.ini"
    #start_ra = '0h0m0s'
    #start_dec = '90d0m0s'
    #single_blra = '8h0m12s'
    #single_bldec = '73d55m20s'
    rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)
    # 创建配置文件对象
    config = configparser.ConfigParser()
    # sdr_freq,sdr_bw,sdr_g,sdr_o,sdr_fs,sdr_fn1,sdr_fn2,sdr_path
    t_bl = threading.Thread( target=run_sdr,
    args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1bl,fn2bl,path,) )

    #################
    # record baseline
    #################
    rc.RaDec_goto(single_blra,single_bldec)
    sleeptime = radec_dist(start_ra,start_dec,single_blra,single_bldec,1)
    time.sleep(sleeptime)

    temp_date = datetime.now()
    obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
    t_bl.start()
    rc.tracking(single_blra,single_bldec,obstime/3600)
    t_bl.join()


    # 添加一个节和选项
    config['SDR1-bl_'+ Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn1bl,
    'observation_mode': 'single_point', 'observation_time': obser_time, 
    'coord_format': coord_format, 'coords1': single_blra, 'coords2': single_bldec
    }

    config['SDR2-bl_' + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'baseline_filename':fn2bl,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_blra, 'coords2': single_bldec
    }

    # 将配置写入到文件中
    with open(log_name, 'a') as configfile:
        config.write(configfile)
    rc.close_comm()
    
    return True


def obs_data(cfreq,bandwidth,gain,obstime,fftsize,path,fn1,fn2,start_ra,start_dec,single_ra,single_dec):
    temp_date = datetime.now()
    Date = temp_date.strftime('%Y-%m-%d')
    Time = temp_date.strftime('%H:%M:%S')
    #cfreq = 1420e6
    #bandwidth = 1.8e6
    #gain = 50
    #obstime = 300
    #fftsize = 4096
    #path = './data/single/'
    #fn1 = 'SDR1-data.txt'
    #fn2 = 'SDR2-data.txt'
    coord_format = 'EQ: ra_dec'
    log_name = Date + "_SDR_single_point.ini"
    #start_ra = '0h0m0s'
    #start_dec = '90d0m0s'
    #single_ra = '7h34m35s'
    #single_dec = '-27d01m18s'

    rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)

    #################
    # data record
    #################
    t_data = threading.Thread( target=run_sdr,
    args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )

    config = configparser.ConfigParser()
    # record data
    rc.RaDec_goto(single_ra,single_dec)
    sleeptime = radec_dist(start_ra,start_dec,single_ra,single_dec,1)
    time.sleep(sleeptime)

    temp_date = datetime.now()
    obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
    t_data.start()
    rc.tracking(single_ra,single_dec,obstime/3600)
    t_data.join()

    # 添加一个节和选项
    config['SDR1-data_' + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
    }

    config['SDR2-data_' + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_ra, 'coords2': single_dec
    }

    # 将配置写入到文件中
    with open(log_name, 'a') as configfile:
        config.write(configfile)

    rc.close_comm()
    # single observation end.
    return True

def grid_scan(cfreq,bandwidth,npoints,gapdeg,gain,obstime,fftsize,path,fn1,fn2,start_ra,start_dec,center_ra,center_dec):
    # grid observation
    temp_date = datetime.now()
    Date = temp_date.strftime('%Y-%m-%d')
    Time = temp_date.strftime('%H:%M:%S')
    #cfreq = 1420e6
    #bandwidth = 1.8e6
    #npoints=5
    #gapdeg=7
    #gain = 50
    #obstime = 30
    #fftsize = 4096
    #path = './data/grid/'
    #fn1 = 'SDR1-grid_data.txt'
    #fn2 = 'SDR2-grid_data.txt'
    coord_format = 'EQ: ra_dec'
    log_name = Date + "_SDR_grid_points.ini"
    #start_ra = '1h00m54s'
    #start_dec = '7d12m43s'
    #center_ra = '1h7m55s'
    #center_dec = '7d12m52s'

    rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)

    # record data
    rc.RaDec_goto(center_ra,center_dec)
    sleeptime = radec_dist(center_ra,center_dec,start_ra,start_dec,0)
    print(sleeptime)
    #time.sleep(8)
    time.sleep(sleeptime)

    grids_ra,grids_dec = grid_coords(center_ra,center_dec,npoints,gapdeg)


    for i in range(len(grids_ra)):
        config = configparser.ConfigParser()
        temp_date = datetime.now()
        obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')

        rc.RaDec_goto( grids_ra[i],grids_dec[i] )
        time.sleep( int( np.round(gapdeg/1.5)+1 ) )
        t_data = threading.Thread( target=run_sdr,
        args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )
        t_data.start()
        rc.tracking(grids_ra[i],grids_dec[i],obstime/3600)
        t_data.join()

        # 添加一个节和选项
        config['SDR1-data_'+ str(i+1) + "_" + Time ] = {
        'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
        'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
        'observation_mode': 'single_point', 'observation_time': obser_time,
        'coord_format': coord_format, 'coords1': grids_ra[i], 'coords2': grids_dec[i]
        }

        config['SDR2-data_' + str(i+1) + "_" + Time ] = {
        'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
        'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
        'observation_mode': 'single_point', 'observation_time': obser_time,
        'coord_format': coord_format, 'coords1': grids_ra[i], 'coords2': grids_dec[i]
        }

        # 将配置写入到文件中
        with open(log_name, 'a') as configfile:
            config.write(configfile)
    
    rc.close_comm()
    return True

def obs_data_galatic(cfreq,bandwidth,gain,obstime,fftsize,path,fn1,fn2,start_l,start_b,single_l,single_b):
    temp_date = datetime.now()
    Date = temp_date.strftime('%Y-%m-%d')
    Time = temp_date.strftime('%H:%M:%S')
    #cfreq = 1420e6
    #bandwidth = 1.8e6
    #gain = 50
    #obstime = 300
    #fftsize = 4096
    #path = './data/single/'
    #fn1 = 'SDR1-data.txt'
    #fn2 = 'SDR2-data.txt'
    coord_format = 'Galatic: l_b'
    log_name = Date + "_SDR_single_point.ini"
    #start_ra = '0h0m0s'
    #start_dec = '90d0m0s'
    #single_ra = '7h34m35s'
    #single_dec = '-27d01m18s'

    rc = RotatorController(port='/dev/ttyUSB0',daz=0,dalt=90)

    #################
    # data record
    #################
    t_data = threading.Thread( target=run_sdr,
    args=(str(cfreq),str(bandwidth),str(gain),str(obstime),str(fftsize),fn1,fn2,path,) )

    config = configparser.ConfigParser()
    c = SkyCoord(l=start_l, b=start_b, frame=Galactic, unit=(u.deg, u.deg)).icrs
    start_ra = c.ra.to_string(unit=u.hourangle, sep='hms', precision=1, pad=True)
    start_dec = c.dec.to_string(unit=u.degree, sep='dms', precision=1, pad=True)
    c = SkyCoord(l=single_l, b=single_b, frame=Galactic, unit=(u.deg, u.deg)).icrs
    single_ra = c.ra.to_string(unit=u.hourangle, sep='hms', precision=1, pad=True)
    single_dec = c.dec.to_string(unit=u.degree, sep='dms', precision=1, pad=True)
    # record data
    rc.RaDec_goto(single_ra,single_dec)
    sleeptime = radec_dist(start_ra,start_dec,single_ra,single_dec,1)
    time.sleep(sleeptime)

    temp_date = datetime.now()
    obser_time = temp_date.strftime('%Y-%m-%d %H:%M:%S')
    t_data.start()
    rc.tracking(single_ra,single_dec,obstime/3600)
    t_data.join()

    # 添加一个节和选项
    config['SDR1-data_' + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn1,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_l, 'coords2': single_b
    }

    config['SDR2-data_' + Time ] = {
    'center_freq': cfreq, 'sample_rate': bandwidth, 'gain': gain,
    'obslength': obstime, 'fftsize':fftsize, 'data_path': path, 'data_filename':fn2,
    'observation_mode': 'single_point', 'observation_time': obser_time,
    'coord_format': coord_format, 'coords1': single_l, 'coords2': single_b
    }

    # 将配置写入到文件中
    with open(log_name, 'a') as configfile:
        config.write(configfile)

    rc.close_comm()
    # single observation in galatic coordinates end.
    return True
