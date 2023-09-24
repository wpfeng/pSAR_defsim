#! /usr/bin/env python
"""
  opt.py is one more module of pSAR. A series of optimizations will be 
  scheduled in this module.
  
  The first idea is coming for the topographically correlated signals (TCS) 
  in the Chile case.
  
  Significant TCS was tried to handled by a power-law fitting. Orbital error 
  will also be treated together.
  
  Particle Swarm Optimization should also be considered within this module...
  
  Spectral analysis was added on 28/08/2017 by FWP @CCRS/NRCan
  
"""
from scipy.stats import sem, t, norm
from scipy import mean
from scipy.optimize import *
import numpy as np
try:
    import shapely
except:
    print(" WARNING: shapely canot be imported properly!!!")
#
###############################################################################
#
def nt_xsamples(xmin,xmax,binwidth=1):
    #
    # nt, normal distribution
    #
    xsamples = np.arange(xmin,xmax,binwidth)
    binc = [(xsamples[i]+xsamples[i+1])/2 for i in range(len(xsamples)-1)]
    return xsamples,binc
#
###############################################################################
#
def nt(x,mu,sigma,mode='percentage'):
    #
    # the mode does not fully function in the current version
    #
    y = 1/(sigma*np.sqrt(2 * np.pi)) * \
        np.exp( - (x - mu)**2 / (2 * sigma**2))
    #
    return y
#
def ci(data,cl=0.95):
    #
    # return confidence interval (CI) at a certain confidence level(CL)
    #
    n = len(data)
    m = mean(data)
    std_err = sem(data)
    h = std_err * t.ppf((1 + cl) / 2, n - 1)
    #
    start_data = m - h
    end_data   = m + h
    ci = abs(end_data - start_data)
    #
    return ci,start_data,end_data
###############################################################################
#
def line2region(prof,buf=0.2):
    #
    listline = [(prof[i,0],prof[i,1]) for i in range(prof.shape[0])]
    #
    lineshapely = shapely.geometry.LineString(listline)
    bufshapely = lineshapely.buffer(buf)
    #
    x,y = bufshapely.exterior.coords.xy
    return x,y
#
def ptsinpolygon(pts,polygon):
    #
    # polygon ID
    polyshapely = shapely.geometry.Polygon(polygon)
    flagarr = np.ones((pts.shape[0]), dtype=bool)
    #
    for i in range(pts.shape[0]):
        cptshapely  = shapely.geometry.Point(pts[i,0],pts[i,1])
        flag = polyshapely.contains(cptshapely)
        flagarr[i] = flag
    return flagarr
    #
    #
def polyarea(x,y):
    #
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
#
def p2dist(lat1, long1, lat2, long2):
    #
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    radiusEARTH=6373
    degrees_to_radians = np.pi/180.0
 
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
 
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
 
    # Compute spherical distance from spherical coordinates.
 
    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) =
    # sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
 
    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + \
           np.cos(phi1) * np.cos(phi2))
    arc = np.arccos( cos )
 
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    dist_km = arc * radiusEARTH
    return dist_km
#############################################################
def roi_dempowlaw(indem,inphs,func,x0,scale=5,maxfev=3000):
    #
    tmpdem = indem[::scale,::scale]
    tmpphs = inphs[::scale,::scale]
    tmpdem = tmpdem[tmpphs!=0]
    tmpphs = tmpphs[tmpphs!=0]
    tmpphs = tmpphs[tmpdem>0]
    tmpdem = tmpdem[tmpdem>0]
    tmpdem,tmpphs = powlaw_sortdata(tmpdem,tmpphs)
    #
    # print(tmpdem.shape)
    coef,res = est_powlawcoef(tmpdem,tmpphs,func,x0,maxfev=maxfev)
    return coef

#############################################################
def demerror_func(coef,x):
    return coef*x
#############################################################
def demerror_errofunc(coef,x,y):
    return y - demerror_func(coef,x)

#
def demerror2perpb(phs,baseline,maxfev=5000):
    #
    demerro,flag = leastsq(demerror_errofunc,0,args=(baseline,phs),maxfev=maxfev)
    demerro = np.mean(phs / baseline)
    return demerro,phs - baseline * demerro
##############################################################
def powlaw_sortdata(x,y,nump=80):
    #
    outdata = []
    for cx in np.arange(x.min(),x.max(),(x.max()-x.min())/nump):
        x0 = cx - (x.max()-x.min())/nump/2
        x1 = cx + (x.max()-x.min())/nump/2
        cy = y[np.logical_and(x>=x0,x<=x1)]
        if cy.shape[0] > 10:
           outdata.append([cx,np.mean(cy)])
    #
    outdata = np.array(outdata)
    return outdata[:,0],outdata[:,1]
##############################################################
def powlaw_fitfun_curvefit(x,a,b,c):
    return a + b * (x ** c)
#
def powlaw_fitfun_1(coef,x,iskeepcons=False):
    if iskeepcons:
        #
        return coef[1] * (x ** coef[2])
    else:
        return coef[0] + coef[1] * (x ** coef[2])
#
def powlaw_fitfun_2(coef,x):
    return coef[0] + coef[1] * x + coef[2] * (x ** coef[3])
#
def powlaw_errfun_2(coef,x,y):
    return y - powlaw_fitfun_2(coef,x)
#
def powlaw_errfun_1(coef,x,y):
    return y - powlaw_fitfun_1(coef,x)
##############################################################
def est_powlawcoef(x,y,func,x0,maxfev=8000):
    #
    # rawx = np.copy(x)
    # x0 = [y[0],0,1,1]
    coef, flag = leastsq(func,x0,args=(x,y),maxfev=maxfev)
        
    return coef,flag
##############################################################
def dpl(x,c0,c1,c2,iskeepcons=False):
    #
    if iskeepcons:
        return c1 *  x ** c2
    else:
        return c0  + c1 *  x ** c2 
#
def dpl_dem(x,c0,c1,c2,c3):
    #
    return  c0  + c1 * x + c2 * x ** c3 
    #
##############################################################
#
def cal_plcof(x,y,p0=None,func=dpl,isplot=False,maxfev=100000):
    #
    popt, pcov = curve_fit(func, x, y, p0=p0, maxfev=maxfev)
    #  
    return popt,pcov
##############################################################
def hst_phs2dem(phs,dem,samples):
    #
    #
    cdem      = dem[np.where(phs!=0)]
    cphs      = phs[np.where(phs!=0)]
    #
    cdem     = cdem.astype(float)
    hdem, bininfo = np.histogram(cdem, bins=samples)
    dims     = hdem.shape
    hphs     = hdem * 0
    hphs     = hphs.astype(float)
    hphs_std = np.copy(hphs)
    #
    for index in range(dims[0]):
        tmpphs = cphs[np.where((cdem > bininfo[index]) & (cdem < bininfo[index+1]))]
        hphs[index]     = np.average(tmpphs)
        hphs_std[index] = np.std(tmpphs)
        hdem[index]     = np.average(bininfo[index:index+2])
    #
    return hdem,hphs,hphs_std
#
##############################################################
def estplsim(phs,dem,p0=None,func=dpl,maxfev=100000,samples=50):
    #
    hdem,hphs,hphs_std =  hst_phs2dem(phs,dem,samples)    
    #
    popt,pcov = cal_plcof(hdem,hphs,p0=p0,func=func,maxfev=maxfev)
    #
    if len(popt)==3:
       #
       sim = dpl(dem,popt[0],popt[1],popt[2])
    else:
       #
       sim = dpl_dem(dem,popt[0],popt[1],popt[2],popt[3])
    #
    corphs = phs - sim
    corphs[np.where(phs==0)] = 0
    sim[np.where(phs==0)]    = 0
    #
    return corphs,sim
##############################################################

     
