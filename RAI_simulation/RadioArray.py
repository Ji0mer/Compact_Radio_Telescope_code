import numpy as np
import matplotlib.pyplot as plt

class Array:

    def __init__(self,x_array,y_array,localname='Boston',lon=-71,alt=42.215,freq=1.42e9):
        """
        New method:
        x_array is antenna array's x coordinate
        y_array is antenna array's y coordinate        
        u,v,w: uvw coordiunate.
        XYZ: XYZ coordinate.
        Ha_ind: hour angle's hour angle index.
        """
        self.x_array = np.array([x_array]).flatten(order="c")
        self.y_array = np.array([y_array]).flatten(order="c")
        self.xvector = np.array([])
        self.yvector = np.array([])
        self.zvector = np.array([])

        # obervation location name
        self.localname = localname
        # longitude and altitude in degrees. Defult is Boston
        # longitude increase to east.
        self.longitude = lon/180*np.pi
        self.altitude = alt/180*np.pi
        # observation freq
        self.freq = freq
        self.lam = 3e8/self.freq

################################################################
        
    def __str__(self):
        return f"Observation location is ({self.localname}:{np.round(self.longitude*180/np.pi,2)},{np.round(self.altitude*180/np.pi,2)})."

    def beam(self):
        """
        This is for future we have a beam pattern
        """
        return 1
    
    def plot_array(self,save=False):
        """
        Plot the array configure
        At a plane.
        """
        maxx = max( max(self.x_array),abs(min(self.x_array)) )
        maxy = max( max(self.y_array),abs(min(self.y_array)) )
        fig = plt.figure(figsize=(16,16))
        plt.title("Radio Interferometer Configuration",fontsize=30)
        plt.arrow(0,0,0,0.8*maxx,head_length=0.4,head_width=0.1)
        plt.arrow(0,0,0.8*maxy,0,head_length=0.4,head_width=0.1)
        plt.text(0.5,0.82*maxx,"North",fontsize=22)
        plt.text(0.82*maxy,0.5,"East",fontsize=22)
        plt.scatter(self.x_array,self.y_array,s=140,c='g',alpha=0.9)
        plt.xlim(-maxx-0.5,maxx+0.5)
        plt.ylim(-maxy-0.5,maxy+0.5)
        plt.xticks(fontsize=26)
        plt.yticks(fontsize=26)
        plt.xlabel("X position/km",fontsize=26)
        plt.ylabel("Y position/km",fontsize=26)
        plt.grid()
        if save == True:
            plt.savefig("./array.png")
        else:
            pass
        plt.show()

    def array_vector(self):
        """
        get arrays baseline vectors at north_east plane.
        """
        for i in range(len(self.x_array)):
            for j in range(len(self.y_array)):
                if i != j:
                    self.xvector = np.append( self.xvector,[self.x_array[i] - self.x_array[j]] )
                    self.yvector = np.append( self.yvector,[self.y_array[i] - self.y_array[j]] )
                    self.zvector = np.append( self.zvector,0 )
                else:
                    pass
        return self.xvector,self.yvector
    
    
    def plot_vector(self,save=False):
        """
        plot the baseline vector at world space plane

        """
        
        xv,yv = self.array_vector()
        fig = plt.figure(figsize=(16,16))
        plt.title("Baseline Vectors",fontsize=30)
        maxx = max( max(xv),abs(min(xv)) )
        maxy = max( max(yv),abs(min(yv)) )
        plt.arrow(0,0,0,0.8*maxx,head_length=0.4,head_width=0.1)
        plt.arrow(0,0,0.8*maxy,0,head_length=0.4,head_width=0.1)
        plt.text(0.5,0.82*maxx,"North",fontsize=22)
        plt.text(0.82*maxy,0.5,"East",fontsize=22)
        plt.scatter(xv,yv,s=140,c='g',alpha=0.9)
        #plt.xlim(-20,20)
        #plt.ylim(-20,20)
        plt.xticks(fontsize=26)
        plt.yticks(fontsize=26)
        plt.xlabel("X/km",fontsize=26)
        plt.ylabel("Y/km",fontsize=26)
        plt.grid()
        if save == True:
            plt.savefig("./array_points.png")
        else:
            pass
        plt.show(fig)  


    def pos2XYZ(self):
        """
        Trans the array elements into XYZ coords
        """
        self.X = np.array([])
        self.Y = np.array([])
        self.Z = np.array([])
        xyz_x = np.array([0,-np.sin(self.altitude),np.cos(self.altitude)])
        xyz_y = np.array([1,0,0])
        xyz_z = np.array([0,np.cos(self.altitude),np.sin(self.altitude)])
        
        for i in range(len(self.xvector)):
            temp_v = np.array([ self.xvector[i],self.yvector[i],self.zvector[i] ])
            self.X = np.append( self.X,np.dot(temp_v,xyz_x) )
            self.Y = np.append( self.Y,np.dot(temp_v,xyz_y) )
            self.Z = np.append( self.Z,np.dot(temp_v,xyz_z) )

        #return self.X,self.Y,self.Z


    def XYZ2uvw(self,H,delta):
        """
        Because the Earth's curvature is very small so that we can directly use lines instead rad to calculate.
        Notice: The XYZ coordnates is a little similar to hour-angle system.
        H: object Hour angle list
        delta: object dec
        Vector X: H=0, d=0
        Vector Y: H=-6h, d=0
        Vector Z: d = 90

        """
        self.u = np.append( self.u, np.sin(H)*self.X + np.cos(H)*self.Y )
        self.v = np.append( self.v, -np.sin(delta)*np.cos(H)*self.X + np.sin(delta)*np.sin(H)*self.Y + np.cos(delta)*self.Z )
        self.w = np.append( self.w, np.cos(delta)*np.cos(H)*self.X - np.cos(delta)*np.sin(H)*self.Y + np.sin(delta)*self.Z )
        self.Ha_ind = np.append( self.Ha_ind,H*180/(15*np.pi) )

    def transXYZ2uvw(self,delta,Hs=-4,He=4,dt=0.01):
        """
        将阵列里所有的X，Y，Z按照给的时角序列画出轨迹
        """
        self.Ha_ind = np.array([])
        self.u = np.array([])
        self.v = np.array([])
        self.w = np.array([])
    
        delta = delta*np.pi/180
        Hlist = np.arange(Hs,He,dt)*15/180*np.pi
        for j in range(len(Hlist)):
            self.XYZ2uvw(Hlist[j],delta)
        #return self.u,self.v,self.w,self.Ha_ind



        