#%%
import time
from datetime import datetime, timedelta
import numpy as np
import healpy as hp
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
#from astropy.coordinates import get_moon
import matplotlib.pyplot as plt

from pygdsm import GSMObserver2016, GSMObserver

#%%
def during_time(ra,dec,loc,start_time,thresh):
    """
    ra in degree
    dec in degree
    loc is observation location
    start_time in utc
    thresh is a alt angle you set, it will renturn how long
    time the sky object below in this threshold.

    return time unit in hours.
    
    """
    coord = SkyCoord(ra, dec, unit='deg')

    # 计算一天内的时间点，从观测时间开始算起，共计算1440个时间点，每个时间点相差一分钟
    time_range = start_time + np.linspace(0, 1, 1440)*24*60*60*u.second

    # 计算每个时间点的天体地平坐标
    altaz_range = coord.transform_to(AltAz(obstime=time_range, location=loc))

    # 计算每个时间点的高度角
    alt_range = altaz_range.alt

    if len( np.where(alt_range < thresh*u.deg)[0] ) == 0:
        time_diff = 25
    else:
        # 找到第一个高度角低于你设定角度的时间点
        idx = np.where(alt_range < thresh*u.deg)[0][0]

        # 计算时间间隔
        time_diff = (time_range[idx] - obs_time).value*24

    return time_diff



# %%

# 观测参数
timediff = -4
#obs_date = datetime(2023, 4, 3, 12, 30, 4, 441805) 
obs_date = datetime.now()
print(obs_date)
min_alt_deg = 20
beam_size_deg = 10

obs_location = EarthLocation(lon=-71.0589 * u.deg, lat=42.3601 * u.deg, height=100 * u.m)


# 将观测时间转换为UTC
obs_time = Time(obs_date - timedelta(hours=timediff), scale='utc', location=obs_location)

# 创建HEALPix地图和像素分辨率
nside = 2**2
npix = hp.nside2npix(nside)
hpx_map = np.zeros(npix)
resol = hp.pixelfunc.nside2resol(nside, arcmin=True)
print("Pixel resolution in arcmin:", resol/60)

sun_loc = get_sun(obs_time)

radius_deg = 20

# 将太阳坐标转换为SkyCoord对象
sun_coord = SkyCoord(sun_loc.ra.degree, sun_loc.dec.degree, unit='deg', frame='icrs')

# 获取太阳附近10度范围内的像素索引
vec = hp.ang2vec(sun_loc.ra.degree, sun_loc.dec.degree, lonlat=True)
ipix_disc = hp.query_disc(nside, vec, np.radians(radius_deg))

# ra and dec list for the healpix area center coords
ra_list = []
dec_list = []


# 遍历所有HEALPix像素
for ipix in range(npix):
    # 获取像素的赤经和赤纬坐标
    coord = hp.pix2ang(nside, ipix, lonlat=True)
    ra, dec = coord
    
    # 将赤经赤纬坐标转换为地平坐标
    sky_coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
    alt_az_coord = sky_coord.transform_to(AltAz(obstime=obs_time, location=obs_location))

    # 检查高度角是否大于20度，是否在太阳附近
    if alt_az_coord.alt.deg > min_alt_deg:
        if ipix in ipix_disc:
            hpx_map[ipix] = 0
        else:
            hpx_map[ipix] = 1
            ra_list.append(ra)
            dec_list.append(dec)


#使用Mollweide投影绘制HEALPix地图
hp.mollview(hpx_map, title="Observable Sky", coord=['C'], unit='Sky Coverage', rot=180)
plt.show()

#%%

time1= time.time()

sunset_time = during_time(sun_loc.ra.degree,sun_loc.dec.degree,obs_location,obs_time,0)
print(f"Time to sun set is {sunset_time} hours.")

settime_list = []
for i in range( len(ra_list) ):
    temp_set = during_time(ra_list[i],dec_list[i],obs_location,obs_time,min_alt_deg)
    settime_list.append(temp_set)

time2 = time.time()
print(f"spend time is {time2 - time1}.")


#%%

np.savez("./"+obs_date.strftime("%Y-%m-%d_%H:%M:%S")+"_surveyplan.npz",ra = ra_list, dec = dec_list, settime = settime_list)

# %%
# Setup observatory location
(latitude, longitude, elevation) = ('42.3601', '-71.0589', 100)
ov = GSMObserver2016()
ov.lon = longitude
ov.lat = latitude
ov.elev = elevation
ov.date = datetime(2023, 4, 12, 13, 30) - timedelta(hours=timediff)

ov.generate(1420)
d = ov.view(logged=True)

# %%
d = ov.view_observed_gsm()

# %%
