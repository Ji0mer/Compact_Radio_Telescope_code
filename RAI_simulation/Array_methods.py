import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import optimize

def VLA_A():
    xa = np.array([ -0.399073,  -1.31578 ,  -2.64078 ,  -4.32773 ,  -6.3478  ,
        -8.68028 , -11.3093  , -14.2221  , -17.408   ,   0.440906,
         1.44321 ,   2.89167 ,   4.73584 ,   6.94428 ,   9.49414 ,
        12.3683  ,  15.5527  ,  19.0925  ,  -0.035827,  -0.122767,
        -0.257158,  -0.408449,  -0.600083,  -0.821349,  -1.07072 ,
        -1.34704 ,  -1.64926 ])
    ya = np.array([ -0.0303058,  -0.648692 ,  -1.54239  ,  -2.68028  ,  -4.0428   ,
        -5.61613  ,  -7.38908  ,  -9.35376  , -11.5025   ,   0.0358317,
        -0.431519 ,  -1.10686  ,  -1.96672  ,  -2.99641  ,  -4.1851   ,
        -5.5249   ,  -7.00924  ,  -8.50824  ,   0.67504  ,   1.66847  ,
         3.20367  ,   4.93181  ,   7.12044  ,   9.6478   ,  12.4962   ,
        15.652    ,  19.1037   ])
    
    return xa,ya


def Y_arm(Annum1,Annum2,dr,dsr,rotation_angle=-85,arm_angle=120):
    axnum = Annum1
    x_array = [0]
    y_array = [0]
    angle = rotation_angle/180*np.pi
    dang = arm_angle/180*np.pi
    for i in range(axnum):
        if i < int(axnum/2):
            x_array.append((i+1)*dr*np.cos(angle))
            y_array.append(-(i+1)*dr*np.sin(angle))
        else:
            x_array.append( (i-int(axnum/2)+1)*dr*np.cos(-angle+dang) )
            y_array.append( (i-int(axnum/2)+1)*dr*np.sin(-angle+dang) )

            x_array.append( (i-int(axnum/2)+1)*dr*np.cos(-angle+2*dang) )
            y_array.append( (i-int(axnum/2)+1)*dr*np.sin(-angle+2*dang) )

    lowfreqnum = Annum2
    if lowfreqnum == 0:
        pass
    else:
        angle = rotation_angle/180*np.pi
        dang = arm_angle/180*np.pi
        for i in range(lowfreqnum):
            x_array.append( (i+1)*dsr*np.cos(angle) )
            y_array.append( -(i+1)*dsr*np.sin(angle) )
    
            x_array.append( -(i+1)*dsr*np.cos(-angle+dang) )
            y_array.append( (i+1)*dsr*np.sin(-angle+dang) )
    
            x_array.append( -(i+1)*dsr*np.cos(-angle+2*dang) )
            y_array.append( (i+1)*dsr*np.sin(-angle+2*dang) )
            
    return x_array,y_array

def Y_arm_power(Annum1,Annum2,dx,dsx,rotation_angle=-85,arm_angle=120):
    axnum = Annum1
    x_array = [0]
    y_array = [0]
    angle = rotation_angle/180*np.pi
    dang = arm_angle/180*np.pi
    for i in range(axnum):
        if i < int(axnum/2):
            x_array.append( np.exp( (i-1)*dx )*np.cos(angle))
            y_array.append(-np.exp( (i-1)*dx )*np.sin(angle))
        else:
            x_array.append( np.exp( (i-int(axnum/2)-1)*dx )*np.cos(-angle+dang) )
            y_array.append( np.exp( (i-int(axnum/2)-1)*dx )*np.sin(-angle+dang) )

            x_array.append( np.exp( (i-int(axnum/2)-1)*dx )*np.cos(-angle+2*dang) )
            y_array.append( np.exp( (i-int(axnum/2)-1)*dx )*np.sin(-angle+2*dang) )
            
    lowfreqnum = Annum2
    if lowfreqnum == 0:
        pass
    else:
        angle = rotation_angle/180*np.pi
        dang = arm_angle/180*np.pi
        for i in range(lowfreqnum):
            x_array.append((i+1)*dsx*np.cos(angle))
            y_array.append(-(i+1)*dsx*np.sin(angle))

            x_array.append( (i+1)*dsx*np.cos(-angle+dang) )
            y_array.append( (i+1)*dsx*np.sin(-angle+dang) )

            x_array.append( (i+1)*dsx*np.cos(-angle+2*dang) )
            y_array.append( (i+1)*dsx*np.sin(-angle+2*dang) )
            
    return x_array,y_array

def get_CleanBeam(dbm,rate=0.02):
    #clean_beam = np.zeros(dbm.shape,dtype="complex")
    yc = int( len(dbm)/2 )
    
    if yc%2 == 0:
        pass
    else:
        yc = yc +1
    yc1 = int( yc - len(dbm)*rate/2 )
    yc2 = int( yc + len(dbm)*rate/2 )

    xc = int( len(dbm[0])/2 )
    
    if xc%2 == 0:
        pass
    else:
        xc = xc +1
    
    xc1 = int( xc - len(dbm[0])*rate/2 )
    xc2 = int( xc + len(dbm[0])*rate/2 )

    psf = dbm[yc1:yc2+1,xc1:xc2+1]
    
    fig = plt.figure(figsize=(12,12))
    plt.pcolormesh( np.log10(abs(psf)+1) )
    plt.show()
    return psf

def overlapIndices(a1, a2, shiftx, shifty):
    if shiftx >=0:
        a1xbeg=shiftx
        a2xbeg=0
        a1xend=a1.shape[0]
        a2xend=a1.shape[0]-shiftx
    else:
        a1xbeg=0
        a2xbeg=-shiftx
        a1xend=a1.shape[0]+shiftx
        a2xend=a1.shape[0]

    if shifty >=0:
        a1ybeg=shifty
        a2ybeg=0
        a1yend=a1.shape[1]
        a2yend=a1.shape[1]-shifty
    else:
        a1ybeg=0
        a2ybeg=-shifty
        a1yend=a1.shape[1]+shifty
        a2yend=a1.shape[1]

    return ( int(a1xbeg), int(a1xend), int(a1ybeg), int(a1yend) ), ( int(a2xbeg), int(a2xend), int(a2ybeg), int(a2yend) )

        
def Clean_Hogbom(dirty,
           psf,
           window,
           gain,
           thresh,
           niter):
    """
    Hogbom clean

    :param dirty: The dirty image, i.e., the image to be deconvolved

    :param psf: The point spread-function, is the dirty beam

    :param window: Regions where clean components are allowed. If
    True, thank all of the dirty image is assumed to be allowed for
    clean components

    :param gain: The "loop gain", i.e., the fraction of the brightest
    pixel that is removed in each iteration

    :param thresh: Cleaning stops when the maximum of the absolute
    deviation of the residual is less than this value

    :param niter: Maximum number of components to make if the
    threshold "thresh" is not hit
    """
    comps = np.zeros(dirty.shape)
    res=np.array(dirty)
    if window is True:
        window=np.ones(dirty.shape,np.bool_)
    for i in range(niter):
        # max value in dirty map
        mx, my=np.unravel_index(np.fabs(res[window]).argmax(), res.shape)
        # peak value and loop gain
        mval=res[mx, my]*gain
        comps[mx, my]+=mval
        a1o, a2o=overlapIndices(dirty, psf,mx-dirty.shape[0]/2,my-dirty.shape[1]/2)
        #print(a1o,a2o)
        res[a1o[0]:a1o[1],a1o[2]:a1o[3]] = res[a1o[0]:a1o[1],a1o[2]:a1o[3]] - psf[a2o[0]:a2o[1],a2o[2]:a2o[3]]*mval
        sub_dm = dirty - res
        if np.fabs(res).max() < thresh:
            print("reach the thresh.")
            break
    return comps, res, sub_dm


def Clean_HB(dirty,psf,gain,thread,nsteps):
    
    res = dirty
    """
    
    psf = psf/np.max(psf)
    psfx = len(psf[0])
    psfy = len(psf)
    dmx = len(dirty[0])
    dmy = len(dirty)
    blx = int( psfx/2 )
    bux = int( dmx-psfx/2 )
    bly = int( psfy/2 )
    buy = int( dmy-psfy/2 )
    """
    
    ind = 0
    clean_model = np.zeros(res.shape)
    
    while ind < nsteps:
        maxres = np.max(res)
        local = np.where( res == np.max(res) )
        #print(local)
        clean_model[local] = maxres
        res[local] = gain*maxres
        """
        if bly<local[0][0]<buy and blx<local[1][0]<bux:
            resly = local[0][0]-int(psfy/2)
            resuy = local[0][0]+int(psfy/2)
            reslx = local[1][0]-int(psfx/2)
            resux = local[1][0]+int(psfx/2)
            res[resly:resuy,reslx:resux] = res[resly:resuy,reslx:resux] - maxres*gain*psf
        else:
            res[local] = maxres*gain
        #print(maxres)
        """
        ind = ind + 1
    return clean_model,res





def con2d(array1,array2,fv=0,mod="full"):
    output = signal.convolve2d(array1,array2,mode=mod,boundary='fill',fillvalue=fv)
    return output


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-y)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p



#########################################################################################################

    



def dir_cosx(vx,vy,vz):
    """
    calculate direction cos
    if the type = "x", it will return direction cos with x axis.
    if the type = "y", it will return direction cos with y axis.
    if the type = "z", it will return direction cos with z axis.
    """
    r1 = np.sqrt(vx**2+vy**2+vz**2)
    return vx/r1

    
def dir_cosy(vx,vy,vz):
    """
    calculate direction cos
    if the type = "x", it will return direction cos with x axis.
    if the type = "y", it will return direction cos with y axis.
    if the type = "z", it will return direction cos with z axis.
    """
    r1 = np.sqrt(vx**2+vy**2+vz**2)
    return vy/r1
    
def dir_cosz(vx,vy,vz):
    """
    calculate direction cos
    if the type = "x", it will return direction cos with x axis.
    if the type = "y", it will return direction cos with y axis.
    if the type = "z", it will return direction cos with z axis.
    """
    r1 = np.sqrt(vx**2+vy**2+vz**2)
    return vz/r1 
    
def ang2vec(theta,phi,r=1):
    """
    Angle to vector.
    Notice: 
    theta,phi in rad !
    the theta is zenith distance,
    phi is azimuth, east is 0 and increase towards to north.

    """
    z = r*np.cos(theta)
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)

    return x,y,z

def vec2ang(vec_x,vec_y,vec_z):
    """
    向量转化为角度theta,phi.
    """
    # normalized vector
    l = np.sqrt(vec_x**2+vec_y**2+vec_z**2)
    x = vec_x/l
    y = vec_y/l
    z = vec_z/l
    if z == 1:
        theta = 0
        phi = 0
    else:
        theta = np.arccos(z)
        sin_theta = np.sqrt(1-z**2)
        #print("sin theta",sin_theta)
        cos_phi = x/sin_theta
        if cos_phi > 1:
            cos_phi = 1
        #print("cos theta",cos_phi)
        sin_phi = y/sin_theta
        if sin_phi > 1:
            sin_phi = 1
        if sin_phi >=0:
            phi = np.arccos(cos_phi)
        else:
            tem_phi = np.arccos(cos_phi)
            phi = 2*np.pi - tem_phi
        
    return theta,phi


def get_lm_matrix(xarray,yarray,H,dec):
    vu = [1,0,0]
    vv = [0,1,0]
    # calculate the dl and dm.
    dH1 = xarray[2] - H
    ddec1 = yarray[2] - dec
    dH2 = xarray[1] - H
    ddec2 = yarray[1] - dec
    temp_theta,phi1 = vec2ang(dH1,ddec1,0)
    temp_theta,phi2 = vec2ang(dH2,ddec2,0)
    sx1,sy1,sz1 = ang2vec(dH1,phi1)
    sx2,sy2,sz2 = ang2vec(dH2,phi2)
    dl = abs( dir_cosx(sx1,sy1,sz1) - dir_cosx(sx2,sy2,sz2) )
    dm = abs( dir_cosy(sx1,sy1,sz1) - dir_cosy(sx2,sy2,sz2) )
    ds = dl*dm
    
    dHl = xarray - H
    ddecl = yarray - dec
    
    l_matrix = np.zeros((200,200))
    m_matrix = np.zeros((200,200))
    
    for i in range(len(dHl)):
        for j in range(len(ddecl)):
            temp_theta,temp_phi = vec2ang(dHl[i],ddecl[j],0)
            sx,sy,sz = ang2vec(dHl[i],temp_phi)
            l = dir_cosx(sx,sy,sz)
            m = dir_cosy(sx,sy,sz)
            l_matrix[j,i] = l
            m_matrix[j,i] = m
    return l_matrix,m_matrix,ds

"""
def VB_matrix(xarray,yarray,Iarray,H,dec,u,v):
    lm,mm,ds = get_lm_matrix(xarray,yarray,H,dec)
    #complex_matrix = np.zeros( (200,200),dtype='complex' )
    cos_matrix = np.cos( 2*np.pi*(u*lm+v*mm) )
    sin_matrix = np.sin( 2*np.pi*(u*lm+v*mm) )
    complex_matrix = cos_matrix - sin_matrix*1j
    VBm = Iarray*np.sqrt( 1 - lm**2 - mm**2 )*complex_matrix*ds
    VB = sum( sum( abs(VBm) ) )
    return VB
"""


def VB_new(lm,mm,xarray,yarray,Iarray,H,dec,u,v):
    """
    Thie method is to derive for a baseline vector at a moment VB value.

    If you need all uv sampling, plz loop....
    Maybe I will promote it in the future....

    VB = int int An(l,m)I(l,m)/sqrt( 1-l^2-m^2 ) e^( -j 2pi (ul+vm) ) dldm

    Because pixel, we have:

    VB = sum sum An(l,m)I(l,m)/sqrt( 1-l^2-m^2 ) e^( -j 2pi(ul+vm) ) dldm

    """
    dl = lm[0][2] - lm[0][1]
    dm = mm[0][2] - mm[0][1]
    ds = dl*dm

    dHl = xarray - H
    ddecl = yarray - dec
    VB = 0
    for i in range(len(dHl)):
        for j in range(len(ddecl)):
            VB = VB+Iarray[j][i]/np.sqrt( 1-lm[j][i]**2-mm[j][i]**2 )*complex( 
                np.cos( 2*np.pi*(u*lm[j][i]+v*mm[j][i]) ), -np.sin(2*np.pi*(u*lm[j][i]+v*mm[j][i])) )*ds

    return VB

def AI(ul,vl,vbm,l,m,du=0.1,dv=0.1):
    AI = complex(0,0)
    ds = du*dv
    for i in range(len(vl)):
        AI = AI + vbm[i]*complex( np.cos(2*np.pi*(ul[i]*l+vl[i]*m) ), np.sin( 2*np.pi*(ul[i]*l+vl[i]*m) ) )*ds
    return AI


def VB(xarray,yarray,Iarray,H,dec,u,v):
    """
    Thie method is to derive for a baseline vector at a moment VB value.

    If you need all uv sampling, plz loop....
    Maybe I will promote it in the future....

    VB = int int An(l,m)I(l,m)/sqrt( 1-l^2-m^2 ) e^( -j 2pi (ul+vm) ) dldm

    Because pixel, we have:

    VB = sum sum An(l,m)I(l,m)/sqrt( 1-l^2-m^2 ) e^( -j 2pi(ul+vm) ) dldm

    """
    VB = complex(0,0)
    vu = [1,0,0]
    vv = [0,1,0]
    # calculate the dl and dm.
    dH1 = xarray[2] - H
    ddec1 = yarray[2] - dec
    dH2 = xarray[1] - H
    ddec2 = yarray[1] - dec
    temp_theta,phi1 = vec2ang(dH1,ddec1,0)
    temp_theta,phi2 = vec2ang(dH2,ddec2,0)
    sx1,sy1,sz1 = ang2vec(dH1,phi1)
    sx2,sy2,sz2 = ang2vec(dH2,phi2)
    dl = abs( dir_cosx(sx1,sy1,sz1) - dir_cosx(sx2,sy2,sz2) )
    dm = abs( dir_cosy(sx1,sy1,sz1) - dir_cosy(sx2,sy2,sz2) )
    ds = dl*dm

    dHl = xarray - H
    ddecl = yarray - dec
    VB = 0
    for i in range(len(dHl)):
        for j in range(len(ddecl)):
            temp_theta,temp_phi = vec2ang(dHl[i],ddecl[j],0)
            sx,sy,sz = ang2vec(dHl[i],temp_phi)
            l = dir_cosx(sx,sy,sz)
            m = dir_cosy(sx,sy,sz)
            VB = VB+Iarray[j][i]/np.sqrt( 1-l**2-m**2 )*complex( 
                np.cos( 2*np.pi*(u*l+v*m) ), -np.sin(2*np.pi*(u*l+v*m)) )*ds

    return VB


def plot_uvsampling(x,y,z,scale=0.5,alp=1):
    fig = plt.figure(figsize=(14,10))
    plt.scatter(x,y,c=z,s=scale,alpha=alp)
    plt.colorbar()
    plt.show() 
