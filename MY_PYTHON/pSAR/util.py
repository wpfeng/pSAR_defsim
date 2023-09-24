#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 22:16:08 2017

@author: wafeng
#
"""
import subprocess
import os
import datetime
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.path import Path
import pSAR
import pyproj
import scipy.spatial.transform
from scipy.interpolate import griddata
#
def interp_byrsc(inroi,inrsc,method='linear'):
    # 
    # method can be as below
    #　　　　linear
    #        nearest
    #        cubic
    # 
    lonm,latm = pSAR.roipac.rsc_to_llm(inroi+'.rsc')
    data = pSAR.roipac.roi_read(inroi)
    data[np.isnan(data)] = 0
    #
    Olon,Olat = pSAR.roipac.rsc_to_llm(inrsc)
    #
    #
    lon1 = lonm.ravel()
    lat1 = latm.ravel()
    data = data.ravel()
    #
    lon1 = lon1[data!=0]
    lat1 = lat1[data!=0]
    data = data[data!=0]
    ipts = np.vstack([lon1,lat1]).T
    #
    # m,n = Olon.shape[0],Olon.shape[1]
    #
    # x0 = Olon.ravel()
    # y0 = Olat.ravel()
    # opts= np.vstack([x0,y0]).T
    #
    outdata = griddata(ipts,data,(Olon,Olat),method=method,fill_value=np.nan)
    #
    return outdata
#

def sartime2eratime(intime,timefmt,hourflag=None):
    #
    # hourflag , current_hour
    #     0 for current hour
    #     -1 for previous hour
    #     1 for next hour
    #
    dt = pSAR.ts.timestr2datetime(intime,fmt=timefmt)
    minute = dt.minute
    #
    if hourflag is None:
        if minute >= 30:
            hourflag = 1
        else:
            hourflag = 0
        #
    #
    if hourflag == 1:
        dt2 = dt + datetime.timedelta(hours = 1)
    else:
        if hourflag == 0:
           dt2 = dt
        else:
           dt2 = dt + datetime.timedelta(hours = -1)
    #
    #
    outdate = '%04d%02d%02d' % (dt2.year,dt2.month,dt2.day)
    outhour = dt2.hour
    #
    #print(indate,outdate,inhour,outhour)
    #
    return outdate,outhour,minute
#
def enu2llh_m(me,mn,malt,lat_org,lon_org,alt_org,transformer1=None,transformer2=None):
    outllh = np.zeros([me.shape[0],3])
    #
    if transformer1 is None:
       transformer1 = pyproj.Transformer.from_crs(
          {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
          {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
          )
    if transformer2 is None:
        transformer2 = pyproj.Transformer.from_crs(
          {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
          {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
          )
    #
    for i in range(me.shape[0]):
        #
        llh = enu2llh(me[i], mn[i], malt[i], lat_org, lon_org, alt_org,\
                transformer1=transformer1, transformer2=transformer2)
        #
        outllh[i,[1,0,2]] = llh
        #
    #
    return outllh
def llh2enu_m(mlat,mlon,malt,lat_org,lon_org,alt_org):
    #
    transformer = pyproj.Transformer.from_crs(
          {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
          {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
          )
    Menu = np.zeros([np.array(mlat).shape[0],2])
    for i in range(mlat.shape[0]):
        #
        enu = llh2enu(mlat[i], mlon[i], malt[i], lat_org, lon_org, alt_org,transformer=transformer)
        Menu[i,:] = [enu[0],enu[1]]
    #
    return Menu
        # print(enu.shape)
def llh2enu(lat, lon, alt, lat_org, lon_org, alt_org,transformer=None):
    if transformer is None:
       transformer = pyproj.Transformer.from_crs(
          {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
          {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
          )
    #
    x, y, z = transformer.transform( lon,lat,  alt,radians=False)
    #
    x_org, y_org, z_org = transformer.transform(lon_org,lat_org,  alt_org,radians=False)
    vec=np.array([[ x-x_org, y-y_org, z-z_org]]).T
    #
    rot1 =  scipy.spatial.transform.Rotation.from_euler('x', -(90-lat_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  scipy.spatial.transform.Rotation.from_euler('z', -(90+lon_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    #
    rotMatrix = rot1.dot(rot3)
    #
    enu = rotMatrix.dot(vec).T.ravel()
    return enu.T
#
def enu2llh(x,y,z, lat_org, lon_org, alt_org,transformer1=None,transformer2=None):
    #
    if transformer1 is None:
       transformer1 = pyproj.Transformer.from_crs(
          {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
          {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
          )
    if transformer2 is None:
       transformer2 = pyproj.Transformer.from_crs(
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        )
    #
    x_org, y_org, z_org = transformer1.transform( lon_org,lat_org,  alt_org,radians=False)
    ecef_org=np.array([[x_org,y_org,z_org]]).T
    #
    rot1 =  scipy.spatial.transform.Rotation.from_euler('x', -(90-lat_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  scipy.spatial.transform.Rotation.from_euler('z', -(90+lon_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    #
    rotMatrix = rot1.dot(rot3)
    # 
    ecefDelta = rotMatrix.T.dot( np.array([[x,y,z]]).T )
    ecef = ecefDelta+ecef_org
    lon, lat, alt = transformer2.transform( ecef[0,0],ecef[1,0],ecef[2,0],radians=False)
    #
    return [lat,lon,alt]
###############################################################################
def run(cmd):
    output        = subprocess.Popen(cmd,stdout=subprocess.PIPE,\
                                     stderr=subprocess.PIPE,shell=True)
    flag          = output.wait()
    prints,errors = output.communicate()
    #
    prints = pSAR.util.bytestoutf8(prints)
    return prints,flag,errors
#
def fdmin(a,b):
    '''
     python function as used in fortran dmin1()
    ----------
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return np.min(np.array([a,b]))
#
def fdsign(a,b):
    '''
    same fucntion as dsign in fortran
    '''
    if b<0:
        a = a*-1
    return a 
#
def aki_mll2local(mlon,mlat,refll=[0,0],r_earth=6378206.4):
    '''
    a numpy array version 
    ----------
    mlon : TYPE
        DESCRIPTION.
    mlat : TYPE
        DESCRIPTION.
    ref : TYPE, optional
        DESCRIPTION. The default is [0,0].
    r_earth : TYPE, optional
        DESCRIPTION. The default is 6378206.4.

    Returns
    -------
    None.

    '''
    mlon = np.array(mlon,dtype=float)
    mlat = np.array(mlat,dtype=float)
    mx = np.copy(mlon)
    my = np.copy(mlat) 
    for i in range(mlon.shape[0]):
        mx[i],my[i] = aki_ll2local(mlon[i],mlat[i],refll=refll,r_earth=r_earth)
    
    return mx,my

def aki_ll2local(lon1,lat1,refll=[0,0],r_earth=6378206.4):
    '''
    calculation ways learnt from disazi(), a fotran code from pscmp(2019)
    
    '''
    PI = np.pi
    PI2= np.pi * 2.0
    DEGTORAD = PI/180.
    #
    latb=refll[1]*DEGTORAD
    lonb=refll[0]*DEGTORAD
    latc=lat1*DEGTORAD
    lonc=lon1*DEGTORAD
    #
    if lonb<0: 
        lonb = lonb + PI2
    if lonc < 0:
        lonc = lonc + PI2
    #
    b=0.5*PI-latb
    c=0.5*PI-latc
    #
    if lonc > lonb:
        aa=lonc-lonb
        if aa <= PI:
            iangle = 1
        else:
            aa = PI2 - aa
            iangle = -1
        #
    else:
        aa=lonb-lonc
        if aa <= PI:
            iangle = -1
        else:
            aa = PI2 - aa
            iangle = 1
    #
    s=np.cos(b)*np.cos(c)+np.sin(b)*np.sin(c)*np.cos(aa)
    a=np.arccos(fdsign(fdmin(abs(s),1.0),s))
    dis=a*r_earth
    if(a*b*c == 0.0):
          angleb=0.0
          anglec=0.0
    else:
          s=0.50*(a+b+c)
          a=fdmin(a,s)
          b=fdmin(b,s)
          c=fdmin(c,s)
          anglec=2.0*np.arcsin(fdmin(1.0,np.sqrt(np.sin(s-a)*np.sin(s-b)/(np.sin(a)*np.sin(b)))))
          angleb=2.0*np.arcsin(fdmin(1.0,np.sqrt(np.sin(s-a)*np.sin(s-c)/(np.sin(a)*np.sin(c)))))
          if(iangle == 1):
            angleb=PI2-angleb
          else:
            anglec=PI2-anglec
    #
    #
    # cartesian coordinates of the sation
    #
    ynorth=dis*np.cos(anglec)
    xeast=dis*np.sin(anglec)
    return xeast, ynorth
        
    
def ll2local(lon1,lat1,refll=[0,0],rhumb=True,r_earth=6378206.4):
    '''
    convert (lon,lat) to a local coordiantes in km
    r_earth, from clarke 1866 equatorial radius of the earth
    learnt from IDL script
    
    rewritten by FWP, @SYSU, 2020/03/19, in python
    
    '''
    k = np.pi/180.
    lon0 = refll[0]
    lat0 = refll[1]
    #
    if rhumb:
      #
      x1 = (lon1-lon0)*k
      #
      lr0 = lat0 * k
      lr1 = lat1 * k
      #
      # Mercator y coordinates. Avoid computing alog(0).
      y0 = np.log(np.tan(np.pi/4 + lr0 / 2))
      y1 = np.log(np.tan(np.pi/4 + lr1 / 2))
      #
      Az = np.arctan2(x1, y1-y0)
      # s is the angular distance between points, in radians.
      s = (lr1-lr0)/np.cos(Az)
      s[lr1 == lr0] = np.abs(x1[lr1==lr0]) * np.cos(lr0)
      #
      dist =  s * r_earth
      y = dist * np.cos(Az)
      x = dist * np.sin(Az)
    else:
      #
      lat0 = refll[1]
      lon0 = refll[0]
      #
      coslt1 = np.cos(k*lat1)
      sinlt1 = np.sin(k*lat1)
      coslt0 = np.cos(k*lat0)
      sinlt0 = np.sin(k*lat0)
      #
      cosl0l1 = np.cos(k*(lon1-lon0))
      sinl0l1 = np.sin(k*(lon1-lon0))
      #
      # Cos of angle between pnts
      cosc = sinlt0 * sinlt1 + coslt0 * coslt1 * cosl0l1
      # Avoid roundoff problems by clamping cosine range to [-1,1].
      #
      sinc = np.sqrt(1.0 - cosc*cosc)
      #
      cosaz = sinc*0
      tmpV = (coslt0 * sinlt1 - sinlt0*coslt1*cosl0l1)
      cosaz[abs(sinc)>1e-7] = tmpV[abs(sinc)>1e-7] / sinc[abs(sinc)>1e-7]
      sinaz                 = sinc * 0
      sinaz[abs(sinc)>1e-7] = sinl0l1[abs(sinc)>1e-7] * coslt1[abs(sinc)>1e-7]/sinc[abs(sinc)>1e-7]
      cosaz[abs(sinc)<1e-7] = 1.0
      sinaz[abs(sinc)<1e-7] = 0.0
      #
      Az   = np.arctan2(sinaz, cosaz)
      dist = np.arccos(cosc) * r_earth
      #
      #
      y = dist * np.cos(Az)
      x = dist * np.sin(Az)
    #
    return x, y, Az / k


def ll2local_utm(lls,refll,un=None,ul=None):
    #
    if un is None or ul is None:
       #
       if refll is not None:
           x0,y0,un,ul = pSAR.utm_conversion.from_latlon(refll[1],refll[0])
       else:
           _,_,un,ul = pSAR.utm_conversion.from_latlon(lls[:,1].mean(),\
                                                       lls[:,0].mean())
           x0 = 0
           y0 = 0
    else:
        x0 = 0
        y0 = 0
    x,y = pSAR.utm_conversion.lltoutm(lls[:,1],lls[:,0],\
                                      force_zone_number=un,force_zone_letter=ul)
    return x-x0,y-y0,un,ul
#
def tfmtchange(tstr,fmt1,fmt2):
    # function to change format of a time to another format
    # tstr, time str
    #
    dttim = datetime.datetime.strptime(tstr,fmt1)
    return dttim.strftime(fmt2)
#
def loadSIMlist(inlist):
    #
    # inlist should have identical format with the list used in pSIM_dislocSIM.py
    # originally developed by FWP
    #
    datalist = []
    counter  = 0
    with open(inlist,'r') as fid:
        for cline in fid:
            cline = cline.split('\n')[0]
            tmp   = cline.split()
            if (tmp[0] == "#" and tmp[1] == "Mode" and counter == 0):
                counter = 1
            else:
               #
               if len(tmp)<3:
                  print(" %s" % cline)
                  print(" Error: no enough information is given....")
                  return []
               else:
                  mode = tmp[0]
                  pointing = tmp[1]
                  phsfile = tmp[2]
                  if len(tmp)>3:
                      incfile = tmp[3]
                  else:
                      incfile = "NULL"
                  if len(tmp)>4:
                      azifile = tmp[4]
                  else:
                      azifile = "NULL"
                  datalist.append([mode,pointing,phsfile,incfile,azifile])
               #
               counter += 1
    return np.array(datalist)
#
###############################################################################
#
def isnumber(instring):
    #
    try:
        float(instring)
        return True
    except:
        return False
#
def ftp_isfile(infile,ftpsession):
    '''
    check if the input (infile) is file. True for yes and False for no
    
    '''
    currentdir = ftpsession.pwd()
    try:
        ftpsession.cwd(infile)
    except:
        ftpsession.cwd(currentdir)
        return True
    #
    ftpsession.cwd(currentdir)
    return False

def meanstep(indata):
    #
    steps = np.zeros(indata.shape[0]-1)
    for i in range(steps.shape[0]):
        #
        steps[i] = indata[i+1] - indata[i]
    return np.mean(steps)
#
def pts2dist(pts):
    #
    x0,y0,un,ul = pSAR.utm_conversion.from_latlon(np.mean(pts[:,1]),\
                                                  np.mean(pts[:,0]))
    ux,uy = pSAR.utm_conversion.lltoutm(pts[:,1],pts[:,0],\
                                        force_zone_number=un,\
                                        force_zone_letter=ul)
    ux,uy = ux/1000., uy/1000.
    #
    tol_lens = 0.
    for i in range(ux.shape[0]-1):
        p1 = [ux[i],uy[i]]
        p2 = [ux[i+1],uy[i+1]]
        lens = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
        tol_lens += lens
    return tol_lens
#
def pts2azi(pts):
    #
    #
    x0,y0,un,ul = pSAR.utm_conversion.from_latlon(np.mean(pts[:,1]),\
                                                  np.mean(pts[:,0]))
    ux,uy = pSAR.utm_conversion.lltoutm(pts[:,1],pts[:,0],\
                                        force_zone_number=un,\
                                        force_zone_letter=ul)
    ux,uy = ux/1000., uy/1000.
    #
    outs = []
    for i in range(ux.shape[0]-1):
        p1 = [ux[i],uy[i]]
        for j in range(i,ux.shape[0]):
            p2 = [ux[j],uy[j]]
            x,y = p2[0] - p1[0], p2[1] - p1[1]
            ang1 = np.arctan2(x,y) * 180. / np.pi
            lens = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
            outs.append([ang1,lens])
    #
    outs = np.array(outs)
    #
    total_lens = np.sum(outs[:,1])
    # meanazi = -1 * np.sum(outs[:,1]/total_lens * outs[:,0]) + 90.
    meanazi = np.mean(outs[:,0])-90
    #
    if meanazi < 0:
        meanazi = meanazi + 360.
    #
    return meanazi
    
def removelines(in_file,out_file,sr='time'):
    '''
     To remove any line which contains sr, for example "time"
     
    '''
    #tmp_file = in_file+'.updated'
    fid_new = open(out_file,'w')
    with open(in_file,'r') as fid:
        for cline in fid:
            cline = bytestoutf8(cline)
            #
            if sr not in cline[0:4]:
               fid_new.write(cline)
    #
    fid_new.close()
    if os.path.exists(out_file):
       return True
    else:
       return False 
                
def ellips2pts(center,a,b,angle,num=500):
    #
    angs = np.linspace(0,360,num=num)
    x0   = b*np.cos(angs * np.pi / 180.) 
    y0   = a*np.sin(angs * np.pi / 180.) 
    #
    import pSIMP
    rotx,roty = pSIMP.simp_rotax(x0,y0,angle,0,0)
    rotx = rotx + center[0]
    roty = roty + center[1]
    return rotx,roty
#
def pts2ellipse(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)
#
def ellipse(cens,width,height,phi,num=500):
    '''
    Return coordinates of an ellipse based on the equation parameters
    center of ellipse, cens
    width
    length
    phi     rotation angle
    num     number of points to return
    by Wanpeng Feng, @CCRS/NRCan, Canada 2018-02-08
    
    '''
    # x = h + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)  [1]
    # y = k + b*sin(t)*cos(phi) + a*cos(t)*sin(phi)  [2]
    t = np.linspace(0,360,num=num)
    x = cens[0] + width/2  * np.cos(t*np.pi/180.) * np.cos(phi*np.pi/180.) - \
                  height/2 * np.sin(t*np.pi/180.) * np.sin(phi*np.pi/180.)
    y = cens[1] + height/2 * np.sin(t*np.pi/180.) * np.cos(phi*np.pi/180.) + \
                  width/2  * np.cos(t*np.pi/180.) * np.sin(phi*np.pi/180.)
    return x,y
#
def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]
    #
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    #
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)
    #
    center = ellip.center
    width = ellip.width
    length = ellip.height
    angle = ellip.angle
    #
    return center,width,length,angle
#
def ptsINpolygon_Path(poly,pts):
    '''
    
    Parameters
    ----------
    poly : TYPE
        DESCRIPTION.
    pts : TYPE
        DESCRIPTION.

    Returns
    -------
    this is much faster than using shapely with a loop.

    '''
    region = Path(poly)
    return region.contains_points(pts)
#
def ptsINpolygon(polygon,pts):
    #
    from shapely.geometry import Polygon
    from shapely.geometry import Point
    mypoly = Polygon(polygon)
    mflag = []
    for i in range(pts.shape[0]):
        cpt = Point(pts[i,:])
        flag = mypoly.contains(cpt)
        mflag.append(flag)
    return np.array(mflag)

#
def bytestoutf8(t):
    '''
    To convert a bytes-like variable to utf8
    '''
    if isinstance(t,bytes):
        return t.decode("utf-8")
    return t
#
def pi():
    return 299792458.

def log(funcin,logname=None,info='Start'):
    '''
    Log inputs of a script to a local ascii file
    by Wanpeng Feng, @NRcan, 2017-08-21
    '''
    p_func  = funcin[0]
    #
    if logname is None:
       logname = os.path.basename(p_func)
       logname = os.path.splitext(logname)[0]+'.log'
    cnow = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')
    #
    print(" *** %s @%s into %s" % (info,cnow,logname))
    try:
      with open(logname,'a') as fid:
        fid.write('%s: ' % cnow)
        for i,ckey in enumerate(funcin):
           if i < len(funcin)-1:
              fid.write('%s ' % funcin[i])
           else:
              fid.write('%s\n' % funcin[i])
        #
        fid.close()
    except:
        print("ERROR: %s cannot be written in the current folder!!!" % logname)  
    #
    if os.path.exists(logname):
       return True
    else:
       return False
###############################################################################
    

               
    
    
