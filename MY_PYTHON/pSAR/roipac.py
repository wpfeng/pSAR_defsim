#!/usr/bin/env python
"""
  This is my first module to read ROI_PAC style data. All these data should 
  have a header with an extension .rsc.
  
  Developed by Feng, Wanpeng, @NRCan, 2016-03-01
  
"""
import re
import multiprocessing as mp
import os
import sys
from pSAR import ui_rect
import numpy as np
from matplotlib import pylab as plt
#
# Make sure that this can run smoothly in HPC
#
try:
  from shapely.geometry import Polygon
except ImportError:
  1
#
# import scipy.ndimage
#  
from scipy import ndimage 
#
###############################################################################
def datalist2mask(datalist,njob=4):
    #
    results = read_roi_parallel(datalist,njob=njob)
    for i in range(len(datalist)):
        #
        if i == 0:
            mask = results[0]
            mask[np.isnan(mask)] = 0.
            mask[mask!=0] = 1
        else:
            cmask = results[i]
            cmask[np.isnan(cmask)] = 0.
            cmask[cmask!=0] = 1
            mask = mask * cmask
        #
    #
    lonm,latm = rsc_to_llm(datalist[0]+'.rsc')
    return mask,lonm,latm
#
def read_roi_parallel(datalist,njob=4):
    #
    cpool = mp.Pool(processes=njob)
    #
    # below is an example, allowing multiple inputs in multiprocessing
    #
    # ts_pts_partial = partial(ts_pts,pts=pts)
    # cpool = mp.Pool(processes=njob)
    # roi_read
    #
    results = cpool.map(roi_read,datalist)
    cpool.close()
    cpool.join()
    #
    return results
#
def saramp2map(amp,scale=1.,exp=0.35,to255=False):
    #
    if exp == 1:
        outdata = amp * scale
    else:
        outdata = scale * np.power(amp,exp)
    #
    if to255:
        histv = np.histogram(outdata[outdata!=0],bins=50)
        pdf,pdfv = histv[0],histv[1]
        index = np.where(pdf==np.max(pdf))[0]
        vmean = np.mean(pdfv[index])
        scalor = 125.5/vmean
        print(" pSAR: stretching data with a scalor of %f into 0-255" % scalor)
        outdata = outdata * scalor
    #
    outdata[outdata>255.] = 255.
    outdata[outdata<0]    = 0.
    #
    return outdata

def ext2polygon(ext):
    #
    outpoly  = np.zeros([5,2])
    outpoly[0,:] = [ext[0],ext[2]]
    outpoly[1,:] = [ext[0],ext[3]]
    outpoly[2,:] = [ext[1],ext[3]]
    outpoly[3,:] = [ext[1],ext[2]]
    outpoly[4,:] = [ext[0],ext[2]]
    return outpoly
#
def nonzero_ratio(indata):
    #
    dim = indata.shape
    #
    return np.where(indata!=0)[0].shape[0]/(dim[0]*dim[1])
###############################################################################
def quad_matrixsta(indata,nonzero_portion=0.3,model='var'):
    #
    data = indata.ravel()
    data[np.isnan(data)] = 0.
    num_nonzero = np.count_nonzero(data)
    if float(num_nonzero) / data.shape[0] < (1-nonzero_portion):
       return -1
    else:
       if model.upper() == "VAR":
          return np.var(data[np.nonzero(data)])
       else:
          return np.std(data[np.nonzero(data)])
#
###############################################################################
def rscs2steps(in_rscs):
    steps = []
    for rsc in in_rscs:
        info,ext = rsc_read(rsc)
        steps.append(float(info['X_STEP']))
    #
    steps = np.array(steps)
    return steps.min(),steps.max()
#
def rscs2overlap(in_rscs):
    #
    for i in range(len(in_rscs)-1):
        #
        if i == 0:
           flag, tmppoly = rsc2overlap(in_rscs[i],in_rscs[i+1])
        else:
           flag, tmppoly = rsc2overlap(tmppoly,in_rscs[i+1],mode=[1,0])
        #
    #
    return tmppoly
#
def rsc2overlap(in_rsc_1,in_rsc_2,mode=[0,0]):
    #
    # mode provides data types, 0 for rsc and 1 for polygon
    #
    if mode[0] == 0:
       poly1 = rsc2polygon(in_rsc_1)
    else:
       poly1 = in_rsc_1
    if mode[1] == 0:
       poly2 = rsc2polygon(in_rsc_2)
    else:
       poly2 = in_rsc_2
    #
    poly1_shape = Polygon(poly1)
    poly2_shape = Polygon(poly2)
    flag = poly1_shape.intersects(poly2_shape)
    #
    if flag:
       return flag, np.array(\
              poly1_shape.intersection(poly2_shape).exterior.coords.xy).T
    else:
       return flag, None
#
###############################################################################
#
def ginsarlist(in_ginsar_list):
    #
    fid = open(in_ginsar_list,'r')
    rois = []
    rscs = []
    dinfo= []
    for cline in fid:
        cline = cline.split('\n')[0]
        tmp   = cline.split()
        rois.append(tmp[0])
        rscs.append(tmp[0]+'.rsc')
        dinfo.append([float(tmp[1]),float(tmp[2]),float(tmp[3])])
    #
    return np.array(rois), np.array(rscs), np.array(dinfo)
#
################################################################################
#
def list2rois(in_roi_list,rsc=False):
    #
    fid = open(in_roi_list,'r')
    rois = []
    rscs = []
    for croi in fid:
        croi = croi.split('\n')[0]
        #
        rois.append(croi)
        rscs.append(croi+'.rsc')
    # end
    fid.close()
    #
    rois = np.array(rois)
    rscs = np.array(rscs)
    if rsc:
       return rois,rscs
    else:
       return rois
#
###############################################################################    
#
def rscs2commonpoly(rscs):
    #
    # return common areas from all rscs
    # added by Wanpeng Feng, @NRCan, 2017-02-01
    #
    counter = 0
    for crsc in rscs:
      #
      poly = rsc2polygon(crsc)
      #
      minlon, minlat= np.min(poly[:,0]), np.min(poly[:,1])
      maxlon, maxlat= np.max(poly[:,0]), np.max(poly[:,1])
      counter += 1
      if counter == 1:
         #
         outminlon, outminlat = minlon, minlat
         outmaxlon, outmaxlat = maxlon, maxlat
      else:
         if outminlon > minlon:
            outminlon = minlon
         if outminlat > minlat:
            outminlat = minlat
         if outmaxlon < maxlon:
            outmaxlon = maxlon
         if outmaxlat < maxlat:
            outmaxlat = maxlat
      #
    outpoly = np.zeros([5,2])
    outpoly[0,:] = [outminlon,outminlat]
    outpoly[1,:] = [outminlon,outmaxlat]
    outpoly[2,:] = [outmaxlon,outmaxlat]
    outpoly[3,:] = [outmaxlon,outminlat]
    outpoly[4,:] = [outminlon,outminlat]
    #
    return outpoly
#
###############################################################################
#
def rsc2polygon(rsc):
    #
    info, ext = rsc_read(rsc)
    #
    outpoly  = np.zeros([5,2])
    outpoly[0,:] = [ext[0],ext[2]]
    outpoly[1,:] = [ext[0],ext[3]]
    outpoly[2,:] = [ext[1],ext[3]]
    outpoly[3,:] = [ext[1],ext[2]]
    outpoly[4,:] = [ext[0],ext[2]]
    return outpoly
#
###############################################################################
#    
def losvec(inc,azi,mode='range',pointingdir='right'):
    #
    # inc and azi are both in degree
    # 
    inc_rad = inc / 180. * np.pi
    azi_rad = azi / 180. * np.pi
    #
    if pointingdir.upper() == "LEFT":
       vecscalar = -1.
    else:
       vecscalar = 1.
    ###########################################################################
    if mode.upper() == "RANGE":
        n_vec =      np.sin(azi_rad) * np.sin(inc_rad) * vecscalar
        e_vec = -1 * np.cos(azi_rad) * np.sin(inc_rad) * vecscalar
        u_vec =      np.cos(inc_rad)
    else:
        n_vec = np.cos(azi_rad)
        e_vec = np.sin(azi_rad)
        u_vec = n_vec * 0.
    #
    return e_vec,n_vec,u_vec
###############################################################################
#    
def ehdr_read(ehdr,refrsc=None):
    #
    if refrsc is None or not os.path.exists(refrsc):
       ginfo,gext = roipac_info()
    else:
       ginfo,gext = rsc_read(refrsc)
    #
    fid = open(ehdr,'r')
    #
    for cline in fid:
        cline = cline.split('/n')[0]
        tmp   = cline.split()
        # remove the space elements in a list
        tmp   = [x for x in tmp if x]
        ckey  = tmp[0]
        valout= tmp[1]
        # print(" KEY: %s VAL: %s" % (ckey,valout))
        if ckey.upper() == "NROWS":
            ginfo["FILE_LENGTH"] = valout
        if ckey.upper() == "NCOLS":
            ginfo["WIDTH"] = valout
        if ckey.upper() == "ULXMAP":
            ginfo["X_FIRST"] = valout
        if ckey.upper() == "ULYMAP":
            ginfo["Y_FIRST"] = valout
        if ckey.upper() == "XDIM":
            ginfo["X_STEP"] = valout
        if ckey.upper() == "YDIM":
            #ypos = np.float64(valout)
            #ypos = ypos * -1
            ginfo["Y_STEP"] = "-"+valout #("%f" % ypos)
    #
    xmin  = np.float64(ginfo["X_FIRST"])
    xpos  = np.float64(ginfo["X_STEP"])
    width = np.float64(ginfo["WIDTH"]) 
    ypos  = np.float64(ginfo['Y_STEP'])
    xmax  = xmin + (width - 1) * xpos
    ymax  = np.float64(ginfo["Y_FIRST"])
    file_length = np.float64(ginfo["FILE_LENGTH"])
    # ypos  = np.float64(ginfo["Y_STEP"])
    ymin  = ymax + (file_length - 1) * ypos
    #
    ###########################################################################
    ginfo["X_MAX"] = ("%f" % xmax)
    ginfo["Y_MIN"] = ("%f" % ymin)
    #
    return ginfo  
# 
###############################################################################
#
def ehdrinfo():
    #
    # ncols 1375
    # nrows 649
    # cellsize 0.050401
    # xllcorner -130.128639
    # yllcorner 20.166799
    # nodata_value 9999.000000
    # nbits 32
    # pixeltype float
    # byteorder msbfirst
    #
    # add a keyword, sensor by wanpeng Feng,@NRCan, 2017-02-28
    # 
    ehdrinfo = {'NCOLS':     '0',\
                'NROWS':     '0',\
                'XDIM':      '1',\
                'YDIM':      '1',\
                'PIXELTYPE': 'FLOAT',\
                'ULXMAP':    '0',\
                'ULYMAP':    '0',\
                'NBANDS':    '1',\
                'BYTEORDER': 'I',\
                'SENSOR':    'SAR',\
                'HEADING_DEG': '0',\
                'INCIDENCE':   '0',\
                'NBITS':'32',\
                'LAYOUT':    'BIL'}
    #
    return ehdrinfo
###############################################################################
def info_to_ehdrinfo(ginfo,dtype='FLOAT'):
    #
    einfo = ehdrinfo()
    #
    einfo['PIXELTYPE'] = dtype
    #
    if dtype.upper()== 'INT16' or dtype.upper()== 'SIGNEDINT':
        einfo['NBITS'] = '16'
    #
    for ckey in ginfo.keys():
        if ckey.upper() == "WIDTH":
            einfo['NCOLS'] = ginfo['WIDTH']
        if ckey.upper() == "FILE_LENGTH":
            einfo['NROWS'] = ginfo['FILE_LENGTH']
        if ckey.upper() == "X_FIRST":
            einfo['ULXMAP'] = ginfo['X_FIRST']
        if ckey.upper() == "Y_FIRST":
            einfo['ULYMAP'] = ginfo['Y_FIRST']
        if ckey.upper() == "X_STEP":
            einfo['XDIM'] = ginfo['X_STEP']
        if ckey.upper() == "Y_STEP":
            einfo['YDIM'] = str(abs(np.float64(ginfo['Y_STEP'])))
    return einfo
###############################################################################
def einfo_to_ehdrhdr(einfo,outehdr):
    #
    fid = open(outehdr,'w')
    fid.write('%-20s %s\n' % ('BYTEORDER',einfo['BYTEORDER'])) 
    fid.write('%-20s %s\n' % ('LAYOUT',   einfo['LAYOUT']))
    fid.write('%-20s %s\n' % ('NBITS',     einfo['NBITS']))
    fid.write('%-20s %s\n' % ('PIXELTYPE',einfo['PIXELTYPE']))
    fid.write('%-20s %s\n' % ('NCOLS',    einfo['NCOLS']))
    fid.write('%-20s %s\n' % ('NROWS',    einfo['NROWS']))
    fid.write('%-20s %s\n' % ('XDIM',     einfo['XDIM']))
    fid.write('%-20s %s\n' % ('YDIM',     einfo['YDIM']))
    fid.write('%-20s %s\n' % ('ULXMAP',   einfo['ULXMAP']))
    fid.write('%-20s %s\n' % ('ULYMAP',   einfo['ULYMAP']))
    fid.close()
    #
    if os.path.exists(outehdr):
        return True
    else:
        return False
###############################################################################
def envihdr2rsc(inhdr):
    # a function to read ENVI header (.hdr) into pSAR
    # provided by FWP, @SYSU, Zhuhai, 2020/12/22
    #
    info,keys = roipac_info()
    #
    regex_1 = re.compile(r'^(.+?)\s*=\s*({\s*.*?\n*.*?})$',re.M|re.I)
    regex_2 = re.compile(r'^(.+?)\s*=\s*(.*?)$',re.M|re.I)
    #
    #
    with open(inhdr,'r') as fid:
        for line in fid:
            matches=regex_1.findall(line)
            #Remove them from the header
            subhdr=regex_1.sub('',line)
            matches.extend(regex_2.findall(subhdr))
            if (len(matches)>0):
                cinfo = matches[0]
                key,value = cinfo
                # print(key)
                if key=='samples':
                   info['WIDTH'] = int(value)
                if key=='lines':
                   info['FILE_LENGTH'] = int(value)
                if key=='map info':
                   for i,key in enumerate(value.split(',')):
                       if i == 3:
                           info["X_FIRST"] = float(key);
                       if i == 4:
                           info["Y_FIRST"] = float(key);
                       if i == 5:
                           info['X_STEP'] = float(key)
                       if i == 6:
                           info['Y_STEP'] = -1*float(key)
                           #
            #
    #
    outrsc = inhdr.split('.hdr')[0]+'.rsc'
    if not os.path.exists(outrsc):
        info_to_rsc(info,outrsc)
    #
    return info
def roi2geotiff(inroi,outgeotif,dtype='Float32'):
    #
    inphs_rsc= inroi+'.rsc'
    inphs_hdr=inroi.split('.')[0]+'.hdr'
    rsc2ehdr(inphs_rsc,inphs_hdr)
    #
    command_sTr = (' gdal_translate -of GTiff -ot %s %s %s' % (dtype,inroi,outgeotif))
    # print(command_sTr)
    os.system(command_sTr)       
    if os.path.exists(outgeotif):
        return True
    else:
        return False
#
###############################################################################    
def rsc2ehdr(rsc,ehdr,dtype='float'):
    #
    ginfo,ext = rsc_read(rsc)
    einfo = info_to_ehdrinfo(ginfo,dtype=dtype)
    flag  = einfo_to_ehdrhdr(einfo,ehdr)
    return flag
#
###############################################################################
#
def rsc2gext(rsc):
    #
    info,ext = rsc_read(rsc)
    gext = '-R%f/%f/%f/%f' % (ext[0],ext[1],ext[2],ext[3])
    ginc = '-I%s/%s' % (info["X_STEP"],info["Y_STEP"][1::])
    return gext,ginc
###############################################################################
#
def rsc_update(inrsc,refrsc):
    info1,ext1 = rsc_read(inrsc)
    info2,ext2 = rsc_read(refrsc)
    keys = ['X_FIRST','X_STEP','Y_FIRST','Y_STEP','WIDTH','FILE_LENGTH',\
            'X_MAX','Y_MIN','Z_MIN','Z_MAX']
    for i in range(len(keys)):
        info2[keys[i]] = info1[keys[i]]
    #
    info_to_rsc(info2,inrsc)
    return True
#
def ehdr2rsc(ehdr,rsc,refrsc=None):
    #
    # read information from ehdr
    # write info to a rsc file
    #
    ginfo = ehdr_read(ehdr,refrsc=refrsc)
    info_to_rsc(ginfo,rsc)
    #
    if os.path.exists(rsc):
       return True
    else:
       return False         
###############################################################################
#       
def array_2_info(a):
    #
    info,keys = roipac_info()
    info['WIDTH'] = str(a.shape[1])
    info['FILE_LENGTH'] = str(a.shape[0])
    info['Y_FIRST'] = str(a.shape[0]-1)
    return info
#
###############################################################################
#        
def roipac_info():
    #
    # add "sensor" in default by Wanpeng Feng, @NRCan, 2017-02-28
    # add "POINTINGDIR" in default by Wanpneg Feng, @CCRS/NRCan, 2017-04-10
    #
    info = {}
    info['ACQTIME']     = 'NONE'
    info['MTIME']       = 'NONE'
    info['STIME']       = 'NONE'
    info['SUPERMASTER'] = 'NONE'
    info['REFLON']      = 'NONE'
    info['REFLAT']      = 'NONE'
    info['MODE']        = 'NONE'
    info['LOOP_FILES']  = 'NONE'
    info['LOOP_SIGNS']  = 'NONE'
    info['Z_SHIFT']     = 'NONE'
    info['WIDTH']       = '0'
    info['FILE_LENGTH'] = '0'
    info['X_MAX']       = '0'
    info['Y_MIN']       = '0'
    info['X_STEP']      = '1'
    info['Y_STEP']      = '-1'
    info['X_FIRST']     = '0'
    info['Y_FIRST']     = '0'
    info['RNG_SPACING'] = '0'
    info['AZI_SPACING'] = '0'
    info['INCIDENCE']   = '0'
    info['PERPB']       = '0'
    info['TEMPB']       = '0'
    info['WAVELENGTH']  = '0.05554'
    info['HEADING_DEG'] = '0'
    info['SENSOR']      = 'SAR'
    info['SENSOR_M']    = 'SAR'
    info['SENSOR_S']    = 'SAR'
    info['X_UNIT']      = 'DEG'
    info['Z_MIN']       = 'NONE'
    info['Z_MAX']       = 'NONE'
    #
    # A flag showing the unit of the data
    # this can be 'm','cm','mm','rad',or 'degree','pixel'
    #
    info['UNIT']        = 'NONE'
    info['Y_UNIT']      = 'DEG'
    info['MASTER']      = "NONE"
    info['SLAVE']       = "NONE"
    info['POL_MASTER']  = 'NONE'
    info['POL_SLAVE']   = 'NONE'
    info['POLARISATION'] = "NONE"
    info["POINTINGDIR"]  = "NONE"
    info["PROJECTION"]   = "LL"
    info['TRACK'] = 'NONE'
    #
    #
    return info,info.keys()
#
def listslice(inlist,inslice):
    #
    outlist = [inlist[i] for i in inslice]
    return outlist
###############################################################################
def roi_to_xyz(roi,xyzfile,scale=1,of='b',nozero=False,isnan=False,ovrscale=1):
    #
    roi_rsc   = roi+'.rsc'
    lonm,latm = rsc_to_llm(roi_rsc)
    fmt  = roi_to_fmt(roi)
    data = roi_read(roi,dtype=fmt)
    if isnan:
       data[np.isnan(data)] = 0.
    
    #
    # To allow oversampling by giving an integer factor, ovrscale which is greater than 1
    # by Wanpeng Feng, @NRCan, 2017-02-21
    #
    if ovrscale > 1:
       #
       #
       # data = scipy.ndimage.zoom(data,ovrscale,order=1)
       # latm = scipy.ndimage.zoom(latm,ovrscale,order=1)
       # lonm = scipy.ndimage.zoom(lonm,ovrscale,order=1)
       #
       data = ndimage.zoom(data,ovrscale,order=1)
       latm = ndimage.zoom(latm,ovrscale,order=1)
       lonm = ndimage.zoom(lonm,ovrscale,order=1)
    #
    ##
    # Downsampling
    data = data[::scale,::scale]
    latm = latm[::scale,::scale]
    lonm = lonm[::scale,::scale]
    #
    lonm = lonm.flatten()
    latm = latm.flatten()
    data = data.flatten()
    #
    nump   = lonm.shape[0]
    outxyz = np.concatenate((np.reshape(lonm,[nump,1]),\
                             np.reshape(latm,[nump,1]),\
                             np.reshape(data,[nump,1])),axis=1)
    #
    if nozero:
        outxyz = outxyz[np.nonzero(outxyz[:,2]),:]
        outxyz = outxyz[0,:,:]
    #
    if of == "b":
       # output bindary
       outxyz.tofile(xyzfile,sep="")
    else:
       # np.savetxt(xyzfile,outxyz,fmt='%f %f %f')
       fid = open(xyzfile,'w')
       for index in range(outxyz.shape[0]):
           fid.write("%f %f %f\n" % (outxyz[index,0],outxyz[index,1],outxyz[index,2]))
       fid.close()
    #
###############################################################################
#
def roi_math(roi,outroi,muls=1.):
    #
    #
    fmt     = roi_to_fmt(roi)
    roidata = roi_read(roi,dtype=fmt)
    outdata = np.copy(roidata)
    #
    # print(" Factor: %f ROI: %s" % (muls,roi))
    #
    outdata               = outdata * muls
    outdata[roidata == 0] = 0.
    #
    #
    roi_write(outdata, outroi)
    rsc       = roi+'.rsc'
    info, ext = rsc_read(rsc)
    info_to_rsc(info,outroi+'.rsc')
    #
    if (os.path.exists(outroi) and os.path.exists(outroi + '.rsc')):
       return True
    else:
       return False
#       
###############################################################################
#
def roireplaceZERO(roi,zerov,outroi):
    #
    # replace zero by a given value <zerov>
    #
    fmt     = roi_to_fmt(roi)
    roidata = roi_read(roi,dtype=fmt)
    #
    adtype=str(roidata.dtype)
    #
    if 'int' in adtype.lower():
        if np.isnan(zerov):
            zerov = -99999
        else:
            zerov = int(zerov)
    roidata[np.isnan(roidata)] = zerov
    # 
    roidata[roidata==0] = zerov
    roi_write(roidata,outroi)
    rsc = roi+'.rsc'
    info,ext = rsc_read(rsc)
    info_to_rsc(info,outroi+'.rsc')
    #
    return outroi
###############################################################################    
#
def roi_rewrap(roi,rewraproi,minv=-5,maxv=5,scale=1,zerov=None):
    #
    rsc = roi+'.rsc'
    if os.path.exists(rsc) is not True:
       print(" ERROR: %s cannot be found. Check inputs first!!!" % rsc)
       return False
    #
    info,ext = rsc_read(rsc)
    # 
    fmt      = roi_to_fmt(roi)
    roidata  = roi_read(roi,dtype=fmt)
    rewrapdata = wrap(roidata,minv=minv,maxv=maxv,scale=scale)
    #
    # Convert 64-bit float to float32
    # No special declaration, all roi_pac files will be treated as float-32 
    #
    rewrapdata = rewrapdata.astype('float32')
    if zerov is not None:
        rewrapdata[rewrapdata==0] = zerov
        rewrapdata[np.isnan(rewrapdata)] = zerov
    #
    roi_write(rewrapdata,rewraproi)
    #
    rewraproi_rsc = rewraproi+'.rsc'
    info_to_rsc(info,rewraproi_rsc)
    #
    if os.path.exists(rewraproi):
       return True
    else:
       return False
###############################################################################
def cutroi2file(inroi,ext,outroi):
    #
    # RSC is assumed existing in the same folder as inroi
    # 
    rsc = inroi+'.rsc'
    idx1,idy1 = lonlat_to_index(rsc,[ext[0],ext[3]])
    idx2,idy2 = lonlat_to_index(rsc,[ext[1],ext[2]])
    #
    # Error handling, if the given region is out of the image
    #
    info,ext  = rsc_read(rsc) 
    width = np.int(info['WIDTH'])
    file_length = np.int(info['FILE_LENGTH'])
    if idx1 < 0:
       print(" WARNING: ext[0] has been smaller than lower boundry of image")
       print("          Fix it at the boundary")
       idx1 = 0
    if idx1 >= width:
       print(" WARNING: ext[0] has been greater than upper boundry of image")
       print("          Fix it at the boundary")
       idx1=width - 1
       idx2=width - 1
    if idx2 >= width:
       print(" WARNING: ext[1] has been greater than upper boundry of image")
       idx2 = width - 1
    if idy1 < 0:
       print(" WARNING: ext[2] has been smaller than upper boundry of image")
       idy1 = 0
    if idy2 >= file_length:
       print(" WARNING: ext[3] has been greater than upper boundry of image")
       idy2 = file_length - 1
    #
    lon0,lat0 = index_to_lonlat(rsc,[idx1,idy1])
    lon1,lat1 = index_to_lonlat(rsc,[idx2,idy2])
    fmt       = roi_to_fmt(inroi)
    #
    # read raw data
    data      = roi_read(inroi,dtype=fmt)
    cdata     = data[idy1:idy2+1,idx1:idx2+1]
    #
    # Updating info
    cinfo = info.copy()
    cinfo["WIDTH"]       = str(cdata.shape[1])
    cinfo['FILE_LENGTH'] = str(cdata.shape[0])
    cinfo['X_FIRST']     = str(lon0)
    cinfo['Y_FIRST']     = str(lat0)
    #
    roi_write(cdata,outroi)
    outroi_rsc = outroi+'.rsc'
    info_to_rsc(cinfo,outroi_rsc)

###############################################################################
def roi2los(roi,losroi,unit='m'):
    #
    rsc = roi+'.rsc'
    ginfo,gext = rsc_read(rsc)
    wavelength = float(ginfo['WAVELENGTH'])
    #
    fmt  = roi_to_fmt(roi)
    data = roi_read(roi,dtype=fmt)
    #
    if unit.upper()=="M":
       los  = data * -1. * wavelength / (4. * np.pi)
    if unit.upper()=="CM":
       los  = data * -1. * wavelength / (4. * np.pi) * 100
    if unit.upper()=="MM":
       los  = data * -1. * wavelength / (4. * np.pi) * 1000
    #
    roi_write(los,losroi)
    #
    out_rsc = losroi + '.rsc'
    info_to_rsc(ginfo,out_rsc) 
    #
    if os.path.exists(losroi):
       return True
    else:
       return False 
###############################################################################
def rsc_to_info(rsc):
    return rsc_read(rsc)
#
def info4gmt(roiinfo):
    #
    xmin   = np.float64(roiinfo['X_FIRST'])
    width  = int(roiinfo['WIDTH'])
    xpos   = np.float64(roiinfo['X_STEP'])
    xmax   = xmin + (width - 1) * xpos
    #
    ymax   = np.float64(roiinfo['Y_FIRST'])
    length = int(roiinfo['FILE_LENGTH'])
    ypos   = np.float64(roiinfo['Y_STEP'])
    ymin              = ymax + (length - 1) * ypos
    gext = '-R%f/%f/%f/%f' % (xmin,xmax,ymin,ymax)
    ginc = '-I%f/%f' % (xpos,abs(ypos))
    return gext,ginc

def rsc_read(rsc,ll180=False):
    """
       Read rsc into a library variable 
    """
    roiinfo,keys = roipac_info()
    #
    # Add one more default keyword in rsc by Wanpeng Feng, @NRCan, 2017-03-01
    #
    roiinfo["SENSOR"] = 0
    #
    counter_info = {}
    for ckey in roiinfo.keys():
        counter_info[ckey] = 0
    #
    # print(os.path.exists(rsc))
    if not os.path.exists(rsc):
       print(" ERROR: %s doen not exist in this folder..." % rsc)
       return False,False
    #
    counter = 0
    with open(rsc,'r') as fid:
      for line in fid:
          #
          counter += 1
          line_str = line.split('\n')[0]
          line_str = line_str.split()
          #
          if len(line_str) > 1 :
              cvalue = line_str[1]
          else:
              cvalue = ""
          #
          # for some unknown reasons, a key in a given .rsc may not be pre-defined before to be used...
          #
          if len(line_str)>0 and line_str[0] not in roiinfo.keys():
             roiinfo.update([(line_str[0], cvalue)])
             counter_info.update([(line_str[0],0)])
             #
          if len(line_str)>0 and counter_info[line_str[0]] == 0:
              roiinfo[line_str[0]] = cvalue
          else:
             if len(line_str)>0 and len(cvalue) > len(roiinfo[line_str[0]]):
                 roiinfo[line_str[0]] = cvalue
          #
          if len(line_str)>0:
             counter_info[line_str[0]] = counter_info[line_str[0]] + 1         
    # 
    if roiinfo['X_FIRST'] == "":
        roiinfo['X_FIRST'] = "0"
    if roiinfo['X_STEP'] == "":
        roiinfo['X_STEP'] = "1"    
    if roiinfo['Y_FIRST'] == "":
        roiinfo['Y_FIRST'] = "0"
    if roiinfo['Y_STEP'] == "":
        roiinfo['Y_STEP'] = "1"      
    xmin   = np.float64(roiinfo['X_FIRST'])
    width  = int(roiinfo['WIDTH'])
    xpos   = np.float64(roiinfo['X_STEP'])
    xmax   = xmin + (width - 1) * xpos
    #
    ymax   = np.float64(roiinfo['Y_FIRST'])
    length = int(roiinfo['FILE_LENGTH'])
    ypos   = np.float64(roiinfo['Y_STEP'])
    ymin              = ymax + (length - 1) * ypos
    extent            = [xmin,xmax,ymin,ymax]
    # print(xmin,xmax)
    roiinfo['X_MAX'] = str(xmax)
    roiinfo['Y_MIN'] = str(ymin)
    #
    #
    for ckey in roiinfo.keys():
        roiinfo[ckey] = str(roiinfo[ckey]).replace(" ","")
    #
    if ll180:
        #
        if float(roiinfo['X_FIRST'])>180:
            roiinfo['X_FIRST'] = np.float64(roiinfo['X_FIRST'])-360
            #
            xmin   = np.float64(roiinfo['X_FIRST'])
            xmax   = xmin + (width - 1) * xpos
            extent            = [xmin,xmax,ymin,ymax]
        
    return roiinfo,extent
###############################################################################
#
def rsc_to_poly(rsc):
    #
    info,ext  = rsc_read(rsc)
    poly      = np.zeros((5,2))
    poly[0,:] = [ext[0],ext[3]]
    poly[1,:] = [ext[0],ext[2]]
    poly[2,:] = [ext[1],ext[2]]
    poly[3,:] = [ext[1],ext[3]]
    poly[4,:] = [ext[0],ext[3]]
    #
    return poly
###############################################################################
def remove_duplication(dpoints):
    #
    flag = dpoints[:,0] + dpoints[:,1] + dpoints[:,0] * dpoints[:,1]
    dims = dpoints.shape
    #
    flag,index = np.unique(flag,return_index=True)
    #
    if dims[0] == flag.shape[0]:
       return dpoints
    else:
       dims = index.shape
       outpoints = np.zeros((dims[0],3))
       for c_index in index:
           cdp = dpoints[np.where(flag==flag[c_index])]
           outpoints[c_index,0:2] = cdp[0,0:2]
           outpoints[c_index,2]   = np.average(cdp[:,2])
       #
       return outpoints       
###############################################################################
def rpoly_intersec(poly1,poly2):
    #
    # poly1 and poly2 should be numpy arrays of 5 x 2
    # two rectangular polygons
    #
    # This is an old version. We shoiuld turn to an external module, shapely now...
    # noted by Wanpeng Feng, @NRCan, 2017-02-01
    #
    minx1, miny1 = np.min(poly1[:,0]),np.min(poly1[:,1])
    maxx1, maxy1 = np.max(poly1[:,0]),np.max(poly1[:,1])
    #
    minx2, miny2 = np.min(poly2[:,0]),np.min(poly2[:,1])
    maxx2, maxy2 = np.max(poly2[:,0]),np.max(poly2[:,1])
    flag = True
    if minx1 > maxx2 or maxx1 < minx2: 
       flag = False
    if miny1 > maxy2 or maxy1 < miny2:
       flag = False
    if minx2 > maxx1 or maxx2 < minx1:
       flag = False
    if miny2 > maxy1 or maxy2 < miny1:
       flag = False
    #
    polygon = np.zeros((5,2))
    #
    if flag:
       minx = (minx1 > minx2) * minx1 + (minx1 <= minx2 ) * minx2
       miny = (miny1 > miny2) * miny1 + (miny1 <= miny2 ) * miny2
       maxx = (maxx1 > maxx2) * maxx2 + (maxx1 <= maxx2 ) * maxx1
       maxy = (maxy1 > maxy2) * maxy2 + (maxy1 <= maxy2 ) * maxy1 
       polygon[:,0] = [minx,maxx,maxx,minx,minx]
       polygon[:,1] = [miny,miny,maxy,maxy,miny]
    #
    return flag,polygon  
########################################################
def sub_roi_data(data,lonm,latm,poly,scale=1):
    #
    data = data[::scale,::scale]
    lonm = lonm[::scale,::scale]
    latm = latm[::scale,::scale]
    #
    minx,maxx = np.min(poly[:,0]),np.max(poly[:,0])
    miny,maxy = np.min(poly[:,1]),np.max(poly[:,1])
    #
    flagx1 = np.logical_and(lonm>minx,lonm<=maxx)
    flagx2 = np.logical_and(latm>=miny,latm<=maxy)
    outdata  = data[np.where(np.logical_and(flagx1,flagx2))] #  & latM>=miny & latM<=maxy)]
    outlon   = lonm[np.where(np.logical_and(flagx1,flagx2))]
    outlat   = latm[np.where(np.logical_and(flagx1,flagx2))]
    #
    return outlon,outlat,outdata
########################################################
def sub_roi(roi,poly,scale=1):
    #
    # extract part of ROI file by given polygon
    # 
    # 
    rsc = roi + ".rsc"
    #
    if os.path.isfile(rsc):
       lonm,latm = rsc_to_llm(rsc)
    else:
       print(" ***ERROR*** " + rsc + ' was not found.')
       return None
    #
    # Detect format of roi, float32 and Int16 can be supported in current version
    #
    fmt  = roi_to_fmt(roi)
    data = roi_read(roi,dtype=fmt, isplot=False)
    #
    outlon,outlat,outdata = sub_roi_data(data,lonm,latm,poly,scale=scale)
    #
    return outlon,outlat,outdata
#######################################################
def points2grd(points,lonm,latm,method='cubic',\
               fill_value=np.nan,rescale=False):
    #
    from scipy.interpolate import griddata
    zgrd = griddata(points[:,0:2],points[:,2],(lonm,latm),\
                    method=method,fill_value=fill_value,\
                    rescale=rescale)
    return zgrd
#######################################################
def rgriddata(points,ext,exp_pos,fill_value=np.nan,method='cubic'):
    #
    # make dicrete points into a gridded data
    #
    from scipy.interpolate import griddata 
    #
    xlim = np.zeros(2)
    ylim = np.copy(xlim)
    xlim[0],xlim[1] = ext[0],ext[1] 
    ylim[0],ylim[1] = ext[2],ext[3]
    #
    x_width         = np.rint((xlim[1] - xlim[0])/exp_pos) + 1 
    y_length        = np.rint((ylim[1] - ylim[0])/exp_pos) + 1
    x               = np.linspace(xlim[0],xlim[1], x_width)
    y               = np.linspace(ylim[1],ylim[0], y_length)
    xgrd,ygrd       = np.meshgrid(x,y)
    # 
    zgrd = griddata(points[:,0:2],points[:,2],(xgrd,ygrd),method=method,\
                    fill_value=fill_value)
    zgrd = zgrd.astype('float32')  

    return zgrd
 
###############################################################################
def rsc_to_lonlat(rsc):
    #
    info,ext = rsc_read(rsc)
    #
    width = int(info['WIDTH'])
    file_length = int(info['FILE_LENGTH'])
    #
    # Create a vector for lon variation
    #
    min_lon = np.float32(info['X_FIRST'])
    max_lon = np.float32(info['X_STEP']) * (width - 1)  + min_lon
    lonV    = np.linspace(min_lon, max_lon , num=width, dtype='float32')
    #
    # Create a vector for lat variation
    #
    max_lat = np.float32(info['Y_FIRST'])
    min_lat = np.float32(info['Y_STEP']) * (file_length - 1) + max_lat
    latV    = np.linspace(max_lat, min_lat,\
                  num=file_length,dtype='float32')
    #
    return lonV, latV
###############################################################################
def info_to_llm(info,model=0,scale=1):
    #
    width = int(info['WIDTH'])
    file_length = int(info['FILE_LENGTH'])
    #
    # Create a vector for lon variation 
    #
    min_lon = np.float32(info['X_FIRST'])
    max_lon = np.float32(info['X_STEP']) * (width - 1)  + min_lon 
    lonV    = np.linspace(min_lon, max_lon , num=width, dtype='float32')
    #
    # Create a vector for lat variation
    #
    max_lat = np.float32(info['Y_FIRST'])
    min_lat = np.float32(info['Y_STEP']) * (file_length - 1) + max_lat
    latV    = np.flipud(np.linspace(min_lat, max_lat,num=file_length,\
                                    dtype='float32'))
    #
    # Create matrices for lon and lat
    #
    lonM    = np.repeat([lonV], file_length , axis=0) #.astype('float32')
    latM    = np.repeat([latV], width       , axis=0).T
    if scale > 1:
       lonM = lonM[::scale,::scale]
       latM = latM[::scale,::scale]
    #
    # print( min_lon, max_lat)
    if model==0:
       #
       return lonM,latM
       #
    else:
       # 
       # Updated by Wanpeng Feng, @NRCan, 2017-02-13
       # now 1D results are allowed to output
       #
       return lonM.ravel(), latM.ravel() 
#
# Create matrix for lon and lat, which should be in the same size as the 
# by Wanpeng Feng, @NRCan, 2016-03-23
#
def rsc_to_llm(rsc,model=0,scale=1,ll180=False):
    #
    info,ext = rsc_read(rsc,ll180=ll180)
    #
    width = int(info['WIDTH'])
    file_length = int(info['FILE_LENGTH'])
    #
    # Create a vector for lon variation 
    #
    min_lon = np.float32(info['X_FIRST'])
    max_lon = np.float32(info['X_STEP']) * (width - 1)  + min_lon 
    lonV    = np.linspace(min_lon, max_lon , num=width, dtype='float32')
    #
    # Create a vector for lat variation
    #
    max_lat = np.float32(info['Y_FIRST'])
    min_lat = np.float32(info['Y_STEP']) * (file_length - 1) + max_lat
    latV    = np.flipud(np.linspace(min_lat, max_lat,num=file_length,\
                                    dtype='float32'))
    #
    # Create matrices for lon and lat
    #
    lonM    = np.repeat([lonV], file_length , axis=0) #.astype('float32')
    latM    = np.repeat([latV], width       , axis=0).T
    if scale > 1:
       lonM = lonM[::scale,::scale]
       latM = latM[::scale,::scale]
    #
    # print( min_lon, max_lat)
    if model==0:
       #
       return lonM,latM
       #
    else:
       # 
       # Updated by Wanpeng Feng, @NRCan, 2017-02-13
       # now 1D results are allowed to output
       #
       return lonM.ravel(), latM.ravel()   
#
###############################################################################
#
def roi_read(roi_file,dtype=None,isplot=0,dims=None,band=1,\
             scale=1,ingamma=False,iscomplex=False,is_z_shift=False):
    #
    """
      Read rmg data from a ROI_PAC like data binary file. 
      Updated by Wanpeng Feng, @NRCan, 2017-02-21
      Now the data format in GAMMA can be suppored since now (2017-02-21)
      
    """
    roi_fig  = roi_file + '.pdf'
    rsc      = roi_file + '.rsc'
    
    #
    if dtype is None:
       # updated by Wanpeng Feng, @NRCan, 2017-02-21
       #
       if ingamma:
           dtype = ">f4"
       else:
           dtype = roi_to_fmt(roi_file,dims=dims)
    #
    if os.path.exists(rsc):
       info,ext = rsc_read(rsc)   
       nx       = np.int(info['WIDTH'])
       ny       = np.int(info['FILE_LENGTH'])
       z_shift  = info['Z_SHIFT']
       #
    else:
       z_shift = 'NONE'
       if dims is None:
          print(" pSAR.ROIPAC: ERROR -> No dimensions are provided")
          sys.exit(-1)
       else:
          nx = dims[0]
          ny = dims[1]
    #
    #####################################################
    # check lines based on file size and width from RSC
    #
    if (dims is not None and len(dims) > 1):
       lines  = roi_to_lines(roi_file,dtype=dtype,band=band,dims=dims[0])
    else:
       lines  = roi_to_lines(roi_file,dtype=dtype,band=band,dims=None)
    #
    if lines != ny:
       #
       print(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
       print(" + WARNING: %s has %d lines, other than %d !!!" % (roi_file,lines,ny))
       print(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
       ny = int(lines)
    # 
    # A bug was fixed by Wanpeng Feng, @NRCan, 2017-01-23
    #
    data = np.fromfile(roi_file,dtype=dtype,count=int(nx*ny*band),sep="");
    if not iscomplex:
       data = data.reshape([ny,band,nx])
    else:
       data = data.reshape([ny,nx,band]) 
       #
    #
    if isplot:
       #
       if dtype == 'float32' or dtype=='f':
          data[data == 0.0] = np.nan
       plt.imshow(data, extent=ext)
       plt.savefig(roi_fig,dpi=720)
       plt.show()
    #
    if not iscomplex:
      if band==1:
         data = data[::scale,0,::scale]
      else:
         data = data[::scale,1,::scale]   
    else:
         pwr  = data[::scale,::scale,0]
         imag = data[::scale,::scale,1]  
         data = np.arctan2(imag,pwr)
    #
    # updated to allow a shift during reading data
    #
    if is_z_shift and z_shift.upper() != "NONE":
        data[data!=0] = data[data!=0] - float(z_shift)
    return data
##########################################################
def img_to_lines(inimage,width,dtype='float32'):
    #
    # calculate number of lines of input data
    # based on the given datatype and width
    #
    if dtype.upper() == "FLOAT32":
       per_pixel = 4
    elif dtype.upper() in ("INT16",'I4','>i4'):
       per_pixel = 2
    else:
       per_pixel = 8
    fsize = os.path.getsize(inimage)
    return fsize / width / per_pixel
##########################################################
def roi_to_lines(inphs,dtype='float32',band=1,dims=None):
    #
    if dtype.upper() in ("FLOAT32",">F4"):
       per_pixel = 4
    elif dtype.upper() in ("INT16",'I4','>I4'):
       per_pixel = 2
    elif dtype.upper() == 'FLOAT64': 
       per_pixel = 8
    elif dtype.upper() == 'INT8':
       per_pixel = 1
    #
    fsize = os.path.getsize(inphs)
    if dims is None:
       rsc   = inphs+'.rsc'
       #
       info,ext = rsc_read(rsc)
       width    = int(info["WIDTH"])
    else:
       width = dims
    #
    # print(dtype,width,per_pixel,fsize)
    #
    lines    = fsize / width / per_pixel / band
    return lines
##########################################################
# 
def img_to_fmt(inphs,nx,ny):
    fsize    = os.path.getsize(inphs)
    pixels   = nx * ny
    per_pixel= fsize / pixels
    if per_pixel == 4:
       #
       fmt = 'float32'
    elif per_pixel == 8:
       fmt = 'float64'
    else:
       fmt = 'Int16'
    return fmt
#
###############################################################################
#
def roi_to_fmt(inphs,dims=None,band=1):
    """
      Detect formation of input data, float32, float64 or Int16
      #
    """
    rsc = inphs+'.rsc'
    # print(" Test RSC: %s " % rsc)
    fsize       = os.path.getsize(inphs)
    if dims is not None:
        width       = dims[0]
        file_length = dims[1]
    else:
        info,ext    = rsc_read(rsc)
        width       = int(info['WIDTH'])
        file_length = int(info['FILE_LENGTH'])
    pixels     = width * file_length * band
    per_pixel   = int(float(fsize) / pixels)
    expect_size = per_pixel * pixels
    if (expect_size != fsize):
        print(" ++++++++++++++++++++++++++++++++++++++++++++++")
        print(" pSAR.roipac: checking the format of the file:")
        print(" Warning: The file may have been modified at the end!!!")
    #
    if per_pixel == 4: 
       fmt = 'float32'
    elif per_pixel == 8: 
       fmt = 'float64'
    elif per_pixel == 2:
       fmt = 'Int16'
       fmt = 'int16'
    else:
       fmt = 'int8'
    #
    return fmt 
#
###############################################################################
#
def view(data, topext=None, wrap_flag=False, minv=-5, maxv=5,isfile=0,scale=1,\
         cptname='RdYlBu'):
    #
    ext = []
    if isfile != 0 and topext is None:
       rsc = data + '.rsc'
       info, ext = rsc_read(rsc)
       fmt = roi_to_fmt(data)
       # print(fmt)
       data = roi_read(data,dtype=fmt) 
       data = data.astype('float32')
    #
    if topext is not None:
          ext = topext
    rawdata = data 
    #
    # Rewrap phase with a given range
    #
    if wrap_flag:
       data = wrap(data,minv=minv,maxv=maxv)
    #
    data[np.where(data==0)] = np.nan
    #
    if len(ext) > 0:  
       #
       # figid=plt.figure()
       # Updated by Wanpeng Feng, @NRCan, 2017-02-13
       # colormap is allowed to be specified...
       #
       plt.imshow(data[::scale,::scale], extent=ext,cmap=plt.get_cmap(cptname))
    else:
       plt.imshow(data)
    #
    #plt.show(block=False)
    plt.show()
    #
    return rawdata
#
###############################################################################
def loop2cycles(data):
    #
    return np.round(data/(2*np.pi))
#
def phs2cycles(data,minv=-np.pi,maxv=np.pi,scale=1):
    #
    re_wrapped_data = wrap(np.copy(data),minv=minv,maxv=maxv,scale=scale)
    #
    return np.round((data - re_wrapped_data)/(maxv-minv))

#
###############################################################################
#
def wrap(data, minv=-3.14, maxv=3.14, ext=[0,0,0,0], scale=1, isplot=0):
    # rewrap a unwrapped interferogram with a given range
    # Wanpeng Feng, @NRCan, 2016-03-16
    """
     To re-wrap data with a given range, [minv, maxv]
     
    """
    # print(" FWP: %f and %f" % (minv,maxv))
    # 
    np.seterr(invalid='ignore')
    #
    tdata = data[::scale,::scale]
    # 
    # a bug fixed by Wanpeng Feng, @NRCan, 2016-06-22
    # inf in data can result in endless loop
    #
    tdata[np.isinf(tdata)] = 0.
    #
    tdata[np.isnan(tdata)] = 0.
    #
    # danger... test only
    # tdata[np.abs(tdata) > 1000] = 0
    #
    # print(tdata.shape)
    #
    xind,yind = np.where(tdata < minv)
    #
    while len(xind)>0:
       tdata[xind,yind] = tdata[xind,yind] + ( maxv - minv )
       xind,yind = np.where(tdata < minv)
    #
    xind,yind = np.where(tdata > maxv)
    while len(xind)>0:
       tdata[xind,yind] = tdata[xind,yind] - ( maxv - minv )
       xind,yind = np.where(tdata > maxv)
    #
    # tdata[np.isnan(tdata)] = np.nan
    #
    tdata[tdata==0]=np.nan
    if isplot == 1:
       view(tdata,ext=ext)
    #
    return tdata
   
###########################################################
def fix_dem(dem_file,dtype='Int16',refvalue=12000):
    """
      Fix data value 
    """
    dem   = roi_read(dem_file,dtype=dtype)
    #
    # Have no idea why np.abs(dem) doesn't work.
    # Whatever, np.fabs can return absolute values smoothly ...
    # by Feng, Wanpeng, @NRCan, 2016-03-03
    #
    index = np.fabs(dem) > refvalue
    dem[index] = 0
    #
    return dem
#   
###############################################################################
#
def llm2info(lonm,latm):
    #
    info,_ = roipac_info()
    info['WIDTH'] = lonm.shape[1]
    info['FILE_LENGTH'] = lonm.shape[0]
    info['X_FIRST'] = lonm[0,0]
    info['Y_FIRST'] = latm[0,0]
    info['X_STEP'] = lonm[0,1] - lonm[0,0]
    info['Y_STEP'] = latm[1,0] - latm[0,0]
    return info
##############################################################################
def roi_write(data,outfile,fix32=False,dtype=None,outrsc=False):
    """
      Write a Numpy array into a simple binary file.
      Input:
      --------------------------------------------
      <dtype> output data type
      <fix32> flag to force data into a float format
      ---
      Force the output format with "dtype"
      by Wanpeng Feng, @NRCan, 2017-03-28
      
    """
    #
    # add a new keywork fix32 to allow to force the output format 
    # in float32 format...
    #
    if dtype is None:
       dtype = "float32"
    #
    if fix32:
       dtype = 'float32'
    #
    data = data.astype(dtype)
    #
    data.tofile(outfile,sep="")
    if outrsc:
        info = array_2_info(data)
        info_to_rsc(info,outfile+'.rsc')
    return None

###########################################################
def roi2xyzdata(roi,eps=10**-4):
    #
    data1   = roi_read(roi)
    data1[np.isnan(data1)] = 0.
    #
    roirsc1 = roi+'.rsc'
    lonm1,latm1 = rsc_to_llm(roirsc1)
    lonm1 = lonm1[np.nonzero(data1)]
    latm1 = latm1[np.nonzero(data1)]
    data1 = data1[np.nonzero(data1)]
    #
    lonm1 = lonm1[np.abs(data1)>eps]
    latm1 = latm1[np.abs(data1)>eps]
    data1 = data1[np.abs(data1)>eps]
    #
    return np.vstack((lonm1,latm1,data1)).T
#
def rois2xyz(roi1,roi2,outxyz):
    #
    data1   = roi_read(roi1)
    data1[data1==np.nan] = 0
    roirsc1 = roi1+'.rsc'
    lonm1,latm1 = rsc_to_llm(roirsc1)
    #
    mlon  = lonm1[data1!=0]
    mlat  = latm1[data1!=0]
    mdata = data1[data1!=0]
    xyz1  = np.vstack((mlon,mlat,mdata))
    xyz1  = xyz1.transpose()
    #
    data2   = roi_read(roi2)
    data2[data2==np.nan] = 0
    roirsc2 = roi2+'.rsc'
    lonm2,latm2 = rsc_to_llm(roirsc2)
    #
    mlon  = lonm2[data2!=0]
    mlat  = latm2[data2!=0]
    mdata = data2[data2!=0]
    xyz2  = np.vstack((mlon,mlat,mdata))
    xyz2  = xyz2.transpose()
    #
    xyz   = np.vstack((xyz1,xyz2))
    roi_write(xyz,outxyz)
    #
    if os.path.exists(outxyz):
        return True
    else:
        return False
###########################################################        
    
def roi_shift(roi,shiftv,opt=None,adds=None):
    #
    data = roi_read(roi)
    cdata= np.copy(data)
    if adds is None:
      #
      data[data==np.nan] = 0
      cdata = data + shiftv
    else:
      cdata[data>0]  = cdata[data>0]  + shiftv
      cdata[data<=0] = cdata[data<=0] + adds
    if opt is not None:
        lonm1,latm1 = rsc_to_llm(roi+'.rsc')
        simdata = lonm1 * opt[0] + latm1*opt[1] + opt[2]
        cdata   = data - simdata
    cdata[data==0] = 0
    roi_write(cdata,roi)
    #
###########################################################
def roi_update(roi1,roi2,exp_pos=0.001,method='linear',scale=2,valper=0.2):
    #
    #
    rsc1 = roi1 + '.rsc'
    rsc2 = roi2 + '.rsc'
    #
    poly1 = rsc_to_poly(rsc1)
    poly2 = rsc_to_poly(rsc2)
    flag,poly3 = rpoly_intersec(poly1,poly2)
    #
    if flag:
       #
       ext = [np.amin(poly3[:,0]),np.amax(poly3[:,0]),np.amin(poly3[:,1]),\
              np.amax(poly3[:,1])]
       #
       lonM1,latM1,outdata1 = sub_roi(roi1,poly3,scale=scale)
       lonM2,latM2,outdata2 = sub_roi(roi2,poly3,scale=scale)
       #
       lonM1 = lonM1[np.where(outdata1 != 0)]
       lonM2 = lonM2[np.where(outdata2 != 0)]
       latM1 = latM1[np.where(outdata1 != 0)]
       latM2 = latM2[np.where(outdata2 != 0)]
       #
       outdata1 = outdata1[np.where(outdata1 != 0)]
       outdata2 = outdata2[np.where(outdata2 != 0)]
       #
       dims = outdata1.shape;
       points1 = np.zeros((dims[0],3))
       points1[:,0] = lonM1
       points1[:,1] = latM1
       points1[:,2] = outdata1;
       #
       dims = outdata2.shape;
       points2 = np.zeros((dims[0],3))
       points2[:,0] = lonM2
       points2[:,1] = latM2
       points2[:,2] = outdata2;
       #
       # reshape input data into a give coordinate
       #
       sim1   = rgriddata(points1,ext,exp_pos,fill_value=np.nan,method=method)
       sim2   = rgriddata(points2,ext,exp_pos,fill_value=np.nan,method=method)
       #
       subd1  = sim1[np.where(np.logical_and(np.logical_not(np.isnan(sim1)),\
                                             np.logical_not(np.isnan(sim2))))]
       subd2  = sim2[np.where(np.logical_and(np.logical_not(np.isnan(sim1)),\
                                             np.logical_not(np.isnan(sim2))))]
       res    = subd1 - subd2
       tmpres = res - np.mean(res)
       absres = np.absolute(tmpres)
       index  = np.argsort(absres)
       dims   = index.shape
       # 
       # select 20% of best fit points
       # 
       valInd = int(dims[0]*valper)
       res    = res[index[0:valInd]]
       #
       # return a constant shift
       #
       consshift = np.mean(res)
       #
       # update roi2
       #
       data = roi_read(roi2)
       data[np.where(data!=0)] = data[np.where(data!=0)] + consshift
       roi_write(data,roi2)
    else:
       print(" ***ERROR*** " + " no overlap found between " + roi1 + ' and ' + roi2)
###########################################################
def info_to_rsc(rsc_info,output_rsc):
    """
    """
    # 
    tmpinfo,tmpkeys = roipac_info()
    #
    # Sort a list alphabetically
    # Updated by FWP, @NRCan, 2017-04-06
    #
    tmpkeys         = sorted(tmpkeys)
    #
    outrsc_id = open(output_rsc,"w")
    #
    for c_key in tmpkeys: #rsc_info.keys():
        flag = 0
        for infokey in rsc_info.keys():
            if infokey.upper() == c_key:
                flag = flag + 1
        if flag > 0:
          #
          outrsc_id.write('%-20s %s\n' % (c_key,rsc_info[c_key]))
    #      
    for infokey in rsc_info.keys():
        #
        flag = 0
        for c_key in tmpkeys:
            #
            if c_key.upper() == infokey.upper():
                flag = flag + 1
        #
        if flag == 0:
            outrsc_id.write('%-20s %s\n' % (infokey,rsc_info[infokey]))
    # 
    outrsc_id.close()
    # np.savetxt(output_rsc, rsc_info, delimiter=" ", fmt="%15s")
###########################################################
def loc_to_ij(lonv,latv,locs):
    #
    dimx = lonv.shape[0]
    dimy = latv.shape[0]
    # print(dimx,dimy)
    diff_x = lonv - locs[0]
    diff_y = latv - locs[1]
    xindex = np.where(diff_x<=0)
    x1     = xindex[0][-1] 
    x2     = xindex[0][-1]+1
    #
    yindex = np.where(diff_y>=0)
    y1     = yindex[0][-1]
    y2     = yindex[0][-1]+1 
    #
    if (x1 >=0 and x1 < dimx and x2 < dimx and x2>=0 ):
       xind = int(np.float32(x1) + (locs[0]-lonv[x1]) / (lonv[x2]-lonv[x1]))
    else:
       xind = np.nan
    if (y1 >=0 and y1 < dimy and y2 < dimy and y2>=0 ):
       yind = int(np.float32(y1) + (latv[y1]-locs[1]) / (latv[y2]-latv[y1]))
    else:
       yind = np.nan
    return xind, yind
###########################################################
def index_to_lonlat(rsc,index):
    #
    rsc_info,ext = rsc_read(rsc)
    width        = int(rsc_info['WIDTH'])
    file_length  = int(rsc_info['FILE_LENGTH'])
    x_range      = np.linspace(ext[0], ext[1] , num=width)
    y_range      = np.linspace(ext[2], ext[3] , num=file_length)
    index[1]     = file_length - index[1] - 1
    #
    return x_range[index[0]], y_range[index[1]]

###########################################################
def rsc_ll2ind(rsc,ref_lonlat):
    #
    info,ext = rsc_read(rsc)
    width,file_length = int(info['WIDTH']),int(info['FILE_LENGTH'])
    idx,idy = lonlat_to_index(rsc,ref_lonlat)
    if (idx >= 0 and idx <= width and idy >= 0 and idy <= file_length):
       return idx,idy
    else:
       return None,None
###########################################################
def ll_to_index_frominfo(info,ext,lonlat):
    #
    # a bug was fixed by FWP, 2020/02/20
    # 
    id_x = lonlat[:,0]*0
    id_y = lonlat[:,0]*0
    xstep        = np.float64(info['X_STEP'])
    ystep        = np.float64(info['Y_STEP'])
    #
    #
    id_x = np.fix((lonlat[:,0] - ext[0])/xstep+0.5)
    id_y = np.fix((lonlat[:,1] - ext[3])/ystep-0.5) + 1
    #
    return id_x,id_y
#    
###############################################################################
#
def lonlat_to_index(rsc,ref_lonlat,ll180=False,isarr=False):
    #
    rsc_info,ext = rsc_read(rsc,ll180=ll180)
    xstep        = np.float64(rsc_info['X_STEP'])
    ystep        = np.abs(np.float64(rsc_info['Y_STEP']))
    #
    if isarr:
      #
      print(xstep,ext)
      id_x = (ref_lonlat[:,0] - ext[0])/xstep+0.5  # np.argmin(x_range - ref_lonlat[0]))
      id_y = (ext[3] - ref_lonlat[:,1])/ystep-0.5 + 1  
      id_x = id_x.astype('int')
      id_y = id_y.astype('int')
      #
    else:
      id_x = int((ref_lonlat[0] - ext[0])/xstep+0.5)  # np.argmin(x_range - ref_lonlat[0]))
      id_y = int((ext[3] - ref_lonlat[1])/ystep-0.5) + 1  #
    #
    return id_x, id_y
#
###############################################################################
#
def roi_maskbyroi(inroi,refroi):
    #
    data = roi_read(inroi)
    dataref = roi_read(refroi)
    data[dataref==0] = 0.
    data[np.isnan(dataref)] = 0.
    # data = np.nan_to_num(data,copy=True)
    #
    roi_write(data,inroi,fix32=True)
    return True
##
def roi_maskdata(data,rsc,ext):
    #
    sd,idy,idx = roi_subgrid_data(data,rsc,ext)
    data[idy[0]:idy[1],idx[0]:idx[1]] = 0
    return data
    #
###########################################################
def rsc_subgrid(rsc,ext):
    #
    id_x1,id_y1 = lonlat_to_index(rsc,[ext[0],ext[3]])
    id_x2,id_y2 = lonlat_to_index(rsc,[ext[1],ext[2]])
    #
    info,ext = rsc_read(rsc)
    file_length = int(info['FILE_LENGTH'])
    width       = int(info['WIDTH'])
    #
    if id_x1 < 0:
        id_x1 = 0
    if id_x2 >= width:
        id_x2 = width - 1
    if id_y1 < 0:
        id_y1 = 0
    if id_y2 >= file_length:
        id_y2 = file_length - 1
        
    if id_x1 >= 0 and id_x2 <= width:
       if id_y1 >= 0 and id_y2 <= file_length:
          return [id_y1,id_y2],[id_x1,id_x2]
       else:
          return None,None
    else:
       return None,None
###########################################################
def roi_to_ts(indata,datalist,off):
    #
    fmt = roi_to_fmt(indata)
    rsc = indata+'.rsc'
    #
    data = roi_read(indata,dtype=fmt)
    #
    nump   = datalist.shape[0]
    outdata = []
    for ci in np.arange(nump):
        ext = [datalist[ci,0] - off,datalist[ci,0] + off, datalist[ci,1] - off,datalist[ci,1] + off]
        subdata, idy, idx = roi_subgrid_data(data,rsc,ext)
        if subdata is not None:
           vmean = np.mean(subdata[subdata!=0])
           vstd  = np.std(subdata[subdata!=0])
        else:
           vmean = np.nan
           vstd  = np.nan
        outdata.append([datalist[ci,0],datalist[ci,1],vmean,vstd])
     #
    outdata = np.array(outdata)
    return outdata
###########################################################
def roi_subgrid_data(data,rsc,ext):
    #
    id_y,id_x = rsc_subgrid(rsc,ext)
    #
    subdata = None
    if id_x is not None and id_y is not None:
       subdata = data[id_y[0]:id_y[1],id_x[0]:id_x[1]] 
    #
    return subdata,id_y,id_x
############################################################
def roi_subgrid(roi,ext):
    #
    rsc = roi + '.rsc'
    if os.path.exists(rsc):
       # print([ext[0],ext[3]])
       fmt = roi_to_fmt(roi)
       data = roi_read(roi,dtype=fmt)
       return roi_subgrid_data(data,rsc,ext)
       #
    else:
       print(" ERROR: %s doesn't exist." % rsc)
       return None,None,None
############################################################
def ui_roi(roi_file,isview=0):
    #
    # Pick up a rectangle region with your mouse
    #
    obj = ui_rect.draw_rect()
    #
    rawdata = view(roi_file,wrap_flag=True,minv=-10,maxv=10,isfile=1,scale=5)
    #
    m_poly = obj.poly_m
    # print(m_poly.shape)
    dims   = m_poly.shape
    rsc = roi_file + '.rsc'
    #
    for row in range(0,dims[0]):
        ext = m_poly[row,:]
        maxx = ext[1]
        minx = ext[0] 
        maxy = ext[3] 
        miny = ext[2]
        #
        [indx1,indy1] = lonlat_to_index(rsc,[minx,maxy])
        [minx,maxy]   = index_to_lonlat(rsc,[indx1,indy1])
        [indx2,indy2] = lonlat_to_index(rsc,[maxx,miny])
        [maxx,miny]   = index_to_lonlat(rsc,[indx2,indy2])
        ext[0]        = minx
        ext[1]        = maxx
        ext[2]        = miny
        ext[3]        = maxy
        m_poly[row,:] = ext
        #
        subdata = rawdata[indy1:indy2,indx1:indx2]
        #
        # May make a quick figure 
        #
        if isview != 0:
           #
           subdata = rawdata[indy1:indy2,indx1:indx2]
           plt.figure(row)
           plt.imshow(subdata,extent=[minx,maxx,miny,maxy])
           plt.title('Plot(ROI:' + str(row) + ')')
           plt.show(block=False)
        #
        subdata[np.isnan(subdata)] = 0
        phs = subdata[np.where(subdata != 0)]
        if row == 0:
           m_value = np.array([[np.min(phs),np.max(phs),np.average(phs),np.var(phs)]])
        else:
           m_value = np.append(m_value,[[np.min(phs),np.max(phs),np.average(phs),np.var(phs)]])
    #
    return m_value,m_poly
######################################################################
def roi_extractprof(roi_file,data_list,dtype='float32',buf=1):
    #
    # import scipy.ndimage
    #
    roi_rsc      = roi_file + '.rsc'
    roi_data     = roi_read(roi_file, dtype=dtype)
    #
    # Read data list
    #
    if isinstance(data_list,str):
       datalist = np.loadtxt(data_list)
    else:
       datalist = data_list
    dims         = datalist.shape
    output       = []
    lonv,latv    = rsc_to_lonlat(roi_rsc)
    #
    for indx in range(0,dims[0]):
        xloc = datalist[indx,0]
        yloc = datalist[indx,1]
        id_x,id_y = loc_to_ij(lonv,latv,[xloc,yloc])
        if (~np.isnan(id_x) and ~np.isnan(id_y)):                         
           cdata     = roi_data[id_y-buf:id_y+buf,id_x-buf:id_x+buf]
           cdata     = cdata[np.nonzero(cdata)]
           if len(cdata)>0:                  
              outval = np.mean(cdata[np.nonzero(cdata)])
           else:
              outval = np.nan
        else:
           outval    = np.nan
        #
        output.append([xloc,yloc,outval])
    #
    output = np.array(output)
    return output
###############################################################################
# Extract displacements along a profile 
# from a ROI_PAC like binary file
#
def roi_lltovalue(roi_file,reflon,reflat,dtype='float32',buf=1):
    #
    # note that buf should be larger than 0
    # by Wanpeng Feng, @CCRS/CCMEO/NRCan, 2016-12-09
    #
    roi_rsc   = roi_file + ".rsc"
    roi_data  = roi_read(roi_file, dtype=dtype)
    id_x,id_y = lonlat_to_index(roi_rsc,[reflon,reflat])
    outv      = roi_data[id_y-buf:id_y+buf,id_x-buf:id_x+buf]
    #
    outv[np.isnan(outv)] = 0
    outv = outv[np.where(outv!=0)]
    #
    if len(outv) > 0:
       #
       outv = np.mean(outv)
    else:
       outv = np.nan
    #
    return outv
#####################################################################
# Extract displacements along a profile
# from a matrix 
#
def roi_lltovalue_rsc(roi_data,rscfile,reflon,reflat,mode='ll'):
    #
    # reflon and reflat can be ll-like coordiantes, or line and samples can also be 
    # ok for that. The latter will need a mode like "PIX", rather than "LL"
    #
    if mode.upper() == "LL":
       id_x,id_y = lonlat_to_index(rscfile,[reflon,reflat])
    else:
       id_x,id_y = reflon, reflat
    #
    return roi_data[id_y,id_x]
    #
######################################################################
def roi_to_profile(roi_file,data_list,dtype='float32'):
    """
     data_list, an ascii file with a series of data points
     (lon, lat) 
    """
    roi_rsc      = roi_file + '.rsc'
    roi_data     = roi_read(roi_file, dtype=dtype)
    #
    # Read data list 
    #
    if isinstance(data_list,str):
       datalist = np.loadtxt(data_list)
    else:
       datalist = data_list
    dims         = datalist.shape
    output       = []
    roidims      = roi_data.shape
    #
    for indx in range(0,dims[0]):
        #
        xloc = datalist[indx,0]
        yloc = datalist[indx,1]
        id_x,id_y = lonlat_to_index(roi_rsc,[xloc,yloc])
        #
        if (id_x >=0 and id_x<roidims[1] and id_y>=0 and id_y< roidims[0]):
            output.append([xloc,yloc,roi_data[id_y, id_x]])
        else:
            output.append([xloc,yloc,np.nan])
            
    #
    output = np.array(output)
    return output 
