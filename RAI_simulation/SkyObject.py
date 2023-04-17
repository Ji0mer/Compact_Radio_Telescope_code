import numpy as np
import matplotlib.pyplot as plt

from Array_methods import vec2ang
from Array_methods import ang2vec
from Array_methods import dir_cosx
from Array_methods import dir_cosy
from Array_methods import dir_cosz


class SimulateSkyObject:
    """
    This is a simulated sky radio source.
    
    """
    def __init__(self,Dec,Ha=0,hsize=0.7,dsize=0.7,hpixnum=513,dpixnum=513,freq=1.42e9):
        """
        Ha is radio source reference points hour angle in hours.
        Dec is radio source reference points Decinaltion in degrees.
        hsize is Ha length in degree.
        dsize is Dec length in degree.
        
        """
        self.Ha = Ha*15/180*np.pi
        self.Dec = Dec/180*np.pi
        # hsize in degree
        self.hsize = hsize
        # dsize in degree
        self.dsize = dsize
        # pixel number
        self.hpn = hpixnum
        self.dpn = dpixnum
        # square pixel 2 angle resolution unit in rad
        self.xpixang = np.pi*(self.hsize/self.hpn)/180
        # square pixel 2 angle resolution unit in rad
        self.ypixang = np.pi*(self.dsize/self.dpn)/180
        # the max frequency the fft could derive is 1 over angular resulotion in rad.
        self.xmaxfft = 1/self.xpixang
        self.ymaxfft = 1/self.ypixang
        
        self.HaArray = np.linspace(Ha - hsize/2,Ha + hsize/2,self.hpn)
        self.DecArray = np.linspace(Dec - dsize/2,Dec + dsize/2,self.dpn)
        self.freq = 1.42e9
        self.lam = 3e8/self.freq
        self.rate = 110.16
        self.pixlam = ( self.HaArray[2] - self.HaArray[1] )/60*self.rate

        
    def update_hadec(self,Ha):
        self.Ha = Ha
        self.HaArray = np.linspace(Ha - hsize/2,Ha + hsize/2,self.hpn)
        self.DecArray = np.linspace(Dec - dsize/2,Dec + dsize/2,self.dpn)
        
    def simulation(self,numlow=3,numup=10,background=2):
        xarray = np.linspace( -self.hsize/2, self.hsize/2,self.hpn)
        yarray = np.linspace( -self.dsize/2, self.dsize/2,self.dpn)
        sourcenum = int( np.random.uniform(numlow,numup) )
        X,Y = np.meshgrid(xarray,yarray)
        Z = np.random.random_sample((self.hpn,self.dpn))*background + 0.1
        for i in range(sourcenum):
            sigma1 = np.random.uniform(0.1*self.hsize,0.3*self.hsize)
            sigma2 = np.random.uniform(0.1*self.hsize,0.3*self.hsize)
            xcc = np.random.uniform(-0.3*self.hsize,0.3*self.hsize)
            ycc = np.random.uniform(-0.3*self.dsize,0.3*self.dsize)
            rho = np.random.uniform(-1,1)
            Z = Z + ( 2*np.pi*sigma1*sigma2*np.sqrt(1 - rho**2) )**(-1)*np.exp( -1/(2-2*rho**2)*( ((X-xcc)**2/sigma1**2)-((2*rho*(X-xcc)*(Y-ycc))/(sigma1*sigma2))+((Y-ycc)**2/sigma2**2) ) )
        return X,Y,Z
    
    def plot_simulation(self,logged=False,save=False):
        """
        Do the simulation.

        """
        X,Y,Z = self.simulation()
        
        if logged == False:
            fig = plt.figure(figsize=(18,16))
            plt.xlabel("Hour Angle (deg)",fontsize=26)
            plt.ylabel("Declination (deg)",fontsize=26)
            plt.xticks(fontsize=24,rotation=45)
            plt.yticks(fontsize=24)
            plt.pcolormesh(X,Y,Z,shading='auto')
            plt.colorbar()
            
        elif logged == "Same":
            fig = plt.figure(figsize=(34,12))
            plt.subplot(1,2,1)
            plt.xlabel("Hour Angle (deg)",fontsize=26)
            plt.ylabel("Declination (deg)",fontsize=26)
            plt.xticks(fontsize=24,rotation=45)
            plt.yticks(fontsize=24)
            plt.pcolormesh(X,Y,Z,shading='auto')
            plt.colorbar()
            
            plt.subplot(1,2,2)
            plt.xlabel("Hour Angle (deg)",fontsize=26)
            plt.ylabel("Declination (deg)",fontsize=26)
            plt.xticks(fontsize=24,rotation=45)
            plt.yticks(fontsize=24)
            plt.pcolormesh(X,Y,np.log(Z+1),shading='auto')
            plt.colorbar()
            
        else:
            fig = plt.figure(figsize=(18,16))
            plt.xlabel("Hour Angle (deg)",fontsize=26)
            plt.ylabel("Declination (deg)",fontsize=26)
            plt.xticks(fontsize=24,rotation=45)
            plt.yticks(fontsize=24)
            plt.pcolormesh(X,Y,np.log(Z+1),shading='auto')
            plt.colorbar()
        
        plt.show(fig)
        
        if save == True:
            return X,Y,Z
        else:
            pass
        
    def plot_FFT(self,logged="same",save=False,retxy=True):
        X,Y,Z = self.simulation()
        zfft = np.fft.fft2(Z)
        shft2 = np.fft.fftshift(zfft)
        xfft = np.fft.fftfreq(self.hpn,1/(2*self.xmaxfft))
        for i in range(len(xfft)):
            if xfft[i] >= 0:
                pass
            else:
                break
        xfft = np.append(xfft[i:],xfft[:i])/1000

        yfft = np.fft.fftfreq(self.dpn,1/(2*self.ymaxfft))
        for i in range(len(yfft)):
            if yfft[i] >= 0:
                pass
            else:
                break
        yfft = np.append(yfft[i:],yfft[:i])/1000
        
        
        Xfft,Yfft = np.meshgrid(xfft,yfft)
        
        fig = plt.figure(figsize=(14,12))
        plt.pcolormesh(Xfft,Yfft,np.log10( abs(shft2) ),shading='auto' )
        plt.xticks(fontsize=24,rotation=45)
        plt.xlabel(r"u(k$\lambda$)",fontsize=26)
        plt.yticks(fontsize=24)
        plt.ylabel(r"v(k$\lambda$)",fontsize=26)
        #plt.colorbar()
        plt.show()
        if retxy == True:
            return xfft,yfft,shft2
        else:
            pass




class Skyobject:
    """
    This is sky object.
    You should input a image or 
    generate a 2-D array to simulate a radio source.

    """
    def __init__(self,image,H=14.053,dec=54.35,xpixnum =2624,ypixnum=1952,xang=56,yang=42,freq=1.42e9):
        """
        init the radio source.

        if using default simulation. 
        The hour angle is -2h.
        dec is in deg.

        Else, input the hour angle and decination 
        of the radio source. 

        """
        self.image = image
        self.freq=1.42e9
        self.Ha = H*15
        self.Dec = dec
        self.hpn = xpixnum
        self.dpn = ypixnum
        # x angle, y angel in seconds
        self.xang = xang
        self.yang = yang
        # square pixel 2 angle resolution unit in rad
        self.xpixang = np.pi*(self.xang/self.hpn)/(60*180)
        # square pixel 2 angle resolution unit in rad
        self.ypixang = np.pi*(self.yang/self.dpn)/(60*180)
        self.xmaxfft = 1/self.xpixang
        self.ymaxfft = 1/self.ypixang
        # list
        self.halist = np.linspace( -self.xang/2,self.xang/2,self.hpn)
        self.declist = np.linspace( -self.yang/2,self.yang/2,self.dpn)
        
    def plot_image(self):
        fig = plt.figure(figsize=(16,12))
        plt.pcolormesh(self.image,cmap="gray")
        xpixlist = np.linspace(0,len(self.halist),5)
        ypixlist = np.linspace(0,len(self.declist),5)
        xpixshow = np.round( ( xpixlist - int(len(self.halist)/2) )*( self.halist[2] - self.halist[1] ),2 )
        ypixshow = np.round( ( ypixlist - int(len(self.declist)/2) )*( self.declist[2] - self.declist[1] ),2 )
        plt.xticks(xpixlist,xpixshow,fontsize=24,rotation=45)
        plt.yticks(ypixlist,ypixshow,fontsize=24)
        plt.xlabel("RA offset",fontsize=24)
        plt.ylabel("Dec offset", fontsize=24)
        #plt.colorbar()
        plt.show()
        
    def get_AnI(self,An):
        """
        AnI is intensity * Antenna power pattern
        """
        AnI = self.image*An
        return AnI
    
    def plot_AnI(self,An):
        image = self.get_AnI(An)
        
        fig = plt.figure(figsize=(16,12))
        plt.pcolormesh(self.image,cmap="gray")
        xpixlist = np.linspace(0,len(self.halist),5)
        ypixlist = np.linspace(0,len(self.declist),5)
        xpixshow = np.round( ( xpixlist - int(len(self.halist)/2) )*( self.halist[2] - self.halist[1] ),2 )
        ypixshow = np.round( ( ypixlist - int(len(self.declist)/2) )*( self.declist[2] - self.declist[1] ),2 )
        plt.xticks(xpixlist,xpixshow,fontsize=24,rotation=45)
        plt.yticks(ypixlist,ypixshow,fontsize=24)
        plt.xlabel("RA offset",fontsize=24)
        plt.ylabel("Dec offset", fontsize=24)
        #plt.colorbar()
        plt.show()
        
        
    def plot_FFT(self,logged="same",save=False,retxy=True):
        zfft = np.fft.fft2(self.image)
        shft2 = np.fft.fftshift(zfft)
        xfft = np.fft.fftfreq(self.hpn,1/(2*self.xmaxfft))
        for i in range(len(xfft)):
            if xfft[i] >= 0:
                pass
            else:
                break
        xfft = np.append(xfft[i:],xfft[:i])/1000

        yfft = np.fft.fftfreq(self.dpn,1/(2*self.ymaxfft))
        for i in range(len(yfft)):
            if yfft[i] >= 0:
                pass
            else:
                break
        yfft = np.append(yfft[i:],yfft[:i])/1000
        
        
        Xfft,Yfft = np.meshgrid(xfft,yfft)
        
        fig = plt.figure(figsize=(19,12))
        plt.pcolormesh(Xfft,Yfft,np.log10( abs(shft2) ),shading='auto' )
        plt.xticks(fontsize=24,rotation=45)
        plt.xlabel(r"u(k$\lambda$)",fontsize=26)
        plt.yticks(fontsize=24)
        plt.ylabel(r"v(k$\lambda$)",fontsize=26)
        plt.colorbar()
        plt.show()
        if retxy == True:
            return xfft,yfft,shft2
        else:
            pass
        
        
        
        
        
        
