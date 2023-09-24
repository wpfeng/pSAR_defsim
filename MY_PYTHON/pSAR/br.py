#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 14:37:30 2016
This is a new module under pSAR to read a binary file with a specific pointer.
This will be very helpful to read values at specific pixel(s) from a very large file.
The approximate pixel location from longlat will be determined based on the pixel 
spacing and top-left coordinates. The pixel index can only be integer at the moment.
No interpretation is considered yet if the pixel location is at small portion of a single
pixel.

In default, a single-band float format file is applied.
psize, 4 for single float-32 band
       8 for double float-32 band
       
#
@author: wafeng
wanpeng.feng@hotmail.com

"""

import struct
import pSAR
import numpy as np
import pGAMMA
###############################################################################
def br_geoprof_std(infile,lons,lats,subwid,psize=4,fmt='f'):
    #
    rsc = infile+'.rsc'
    info,ext = pSAR.roipac.rsc_read(rsc)
    width = int(info['WIDTH'])
    flen  = int(info['FILE_LENGTH'])
    fid = open(infile,'rb')
    #
    output = []
    #
    for ni in range(lons.shape[0]):
        i,j  = pSAR.roipac.lonlat_to_index(rsc,[lons[ni],lats[ni]])
        MI = [i-subwid,i+subwid]
        MJ = [j-subwid,j+subwid]
        outdata = br_block_fid(fid,MI,MJ,fwidth=width,flength=flen,\
                               psize=psize,fmt=fmt)
        if outdata is None:
            outdata = np.zeros[3]
        else:
            outdata = outdata[outdata!=0]
        #
        cdata = [lons[ni],lats[ni],np.mean(outdata),np.std(outdata)]
        output.append(cdata)
    fid.close()
    return np.array(output)
#ij
def br_ll2ij(rsc,lons,lats,ll180=False):
    #
    info,ext = pSAR.roipac.rsc_read(rsc,ll180=ll180)
    #
    wid = int(info['WIDTH'])
    flength = int(info['FILE_LENGTH'])
    #
    i0,j0  = pSAR.roipac.lonlat_to_index(rsc,[lons.min(),lats.min()],ll180=ll180)
    i1,j1  = pSAR.roipac.lonlat_to_index(rsc,[lons.max(),lats.min()],ll180=ll180)
    i2,j2  = pSAR.roipac.lonlat_to_index(rsc,[lons.max(),lats.max()],ll180=ll180)
    i3,j3  = pSAR.roipac.lonlat_to_index(rsc,[lons.max(),lats.min()],ll180=ll180)
    #
    dbI = np.array([i0,i1,i2,i3])
    dbJ = np.array([j0,j1,j2,j3])
    #
    #
    if dbI.min()<0:
        minI = 0
    else:
        minI = dbI.min()
    if dbI.max()>=wid:
        maxI = wid-1
    else:
        maxI = dbI.max()
    #
    if dbJ.min()<0:
        minJ = 0
    else:
        minJ = dbJ.min()
    #
    if dbJ.max()>=flength:
        maxJ = flength-1
    else:
        maxJ = dbJ.max()
    #
    #
    MI = [minI,maxI]
    MJ = [minJ,maxJ]
    return MI,MJ
#
def br_geoblock(infile,lons,lats,psize=4,fmt='f'):
    #
    rsc = infile+'.rsc'
    info,ext = pSAR.roipac.rsc_read(rsc)
    #
    wid = int(info['WIDTH'])
    flength = int(info['FILE_LENGTH'])
    #
    MI,MJ = br_ll2ij(rsc,lons,lats)
    #
    if len(np.where(np.array(MI)<0)[0])>0:
        outdata = np.zeros(1) 
    else:
        #
        outdata = br_block(infile,MI,MJ,flength=flength,width=wid,\
                       psize=psize,fmt=fmt)
    return outdata
#
###############################################################################
#
def br_geosubstd(infile,inlon,inlat,subwid,psize=4,fmt='f'):
    #
    # subwid, buffer width
    #
    rsc = infile+'.rsc'
    info,ext = pSAR.roipac.rsc_read(rsc)
    #
    wid = int(info['WIDTH'])
    flength = int(info['FILE_LENGTH'])
    #
    i,j  = pSAR.roipac.lonlat_to_index(rsc,[inlon,inlat])
    MI = [i-subwid,i+subwid]
    MJ = [j-subwid,j+subwid]
    outdata = br_block(infile,MI,MJ,flength=flength,width=wid,\
                       psize=psize,fmt=fmt)
    #
    # updated by Wanpeng Feng, @CCRS, 2018-04-04
    # a standard deviation of the data is returned as well.
    #
    outdata[np.isnan(outdata)] = 0
    #
    return np.mean(outdata[outdata!=0]),np.std(outdata[outdata!=0])
#
###############################################################################
#            
def brij_id(fid,i,j,width,psize=4,fmt='f'):
    #
    shiftdist= (width*j + i) * psize
    # find location of data
    try:
      # print(" Good @ %d and %d" % (i,j))
      fid.seek(shiftdist,0)
      outvalue = struct.unpack(fmt,fid.read(psize))
    except:
      # print(" Warnning: out of boundary @ %d and %d!!!" % (i,j))
      outvalue = [np.nan]
    #
    return outvalue
#
###############################################################################
def bytestoutf8(t,code='utf-8'):
    '''
    To convert a bytes-like variable to utf8
    code could also be latin-1
    '''
    if isinstance(t,bytes):
        return t.decode(code)
    return t    
    #
#
def br_complexarray2fid(data,fid,big=True):
    #
    if big:
       fmt = '>f'
    else:
       fmt = '<f'
    outdata = np.vstack((np.real(data),np.imag(data))).T
    outdata = np.ravel(outdata)
    outdata = outdata.astype(fmt)
    outdata.tofile(fid)
    return True
    
def br_complex2fid(data,fid,fmt='>ff'):
    '''
    fmt  
    '>ff'  big-endian, two float numbers
    '<ff'  little-endian, two float numbers
    '''
    for cd in data:
        fid.write(struct.pack(fmt,np.real(cd),np.imag(cd)))
    return True
def br_complex_fid(fid,cl,width,bigendian=True):
    #
    if bigendian:
        #
        # format for files in big-endian
        #
        fmt_L = '>{}f'.format(width*2)
    else:
        fmt_L = '={}f'.format(width*2)
    #
    # 
    shiftsize = (cl-1) * width * 8
    fid.seek(shiftsize,0)
    #
    # outdata = np.zeros(dims[0],dtype=np.complex64)
    # 
    data = struct.unpack(fmt_L,fid.read(width*8))
    #
    #
    data = np.reshape(np.array(data),(width,2))
    #
    return data[:,0] + 1j * data[:,1]
#
def br_complex(infile,cl,iscexml=None,bigendian=True):
    '''
    br_complex a function to read a line from a binary file
    The file can be in big-endian.
    At the moment, SLC data are generated by GAMMA. A SLC.par or RSLC.par 
    should be available to provide data dimension...
    
    + Input:
    cl   current line
    
    by Wanpeng Feng, @CCRS/NRCan, 2017-08-27 
    #
    '''
    #
    if iscexml is None:
       parfile = infile+'.par'
       dims    = pGAMMA.gamma_slcpar2dim(parfile)
    else:
       import pISCE
       info = pISCE.iscexml(iscexml)
       dims = [int(info['WIDTH']),int(info['FILE_LENGTH'])]
    fid     = open(infile,'rb')
    data    = br_complex_fid(fid,cl,dims[0],bigendian=bigendian)
    fid.close()
    #
    return data
        
    
def br_block_fid(fid,MI,MJ,fwidth=None,flength=None,psize=4,fmt='f'):
    '''
    '''
    #
    # print(MI,MJ)
    if MJ[1] <MJ[0]:
        MJ[1] = MJ[0]
    #
    indata   = np.zeros([MI[1]-MI[0]+1,MJ[1]-MJ[0]+1],dtype=np.float32)
    #
    if fwidth is None or flength is None:
       print(" br_block_fid: <fwidth> and <flength> are mandatory")
       return None
    #
    elif (np.min(MI)>=0 and np.max(MI)<fwidth and \
          np.min(MJ)>=0 and np.max(MJ)<flength):
    
      #
       cnj = 0
       for nj in range(MJ[0],MJ[1]+1,1):
         #
         shiftdist= (fwidth*nj + MI[0]) * psize
         fid.seek(shiftdist,0)
         #
         for ni in range(MI[1]-MI[0]+1):
            #
            outvalue = struct.unpack(fmt,fid.read(psize))
            indata[ni,cnj] = outvalue[0]
            #
         cnj += 1
    else:
         indata[:,:] = np.nan
    return indata
#
###############################################################################
#    
def br_block(infile,MI,MJ,width=None,psize=4,flength=None,fmt='f'):
    #
    if width is None:
      rsc      = infile+'.rsc'
      info,ext = pSAR.roipac.rsc_read(rsc)
      fwidth   = int(info["WIDTH"])
      flength  = int(info['FILE_LENGTH'])
    else:
      fwidth = width
    fid      = open(infile,'rb')
    indata   = br_block_fid(fid,MI,MJ,fwidth=fwidth,flength=flength,\
                            psize=psize,fmt=fmt)
    #
    fid.close()
    #
    return indata 
#    
###############################################################################
def brij_fid(fid,i,j,fwidth=None,psize=4,fmt='f'):
    #
    # print(i,j)
    if fwidth is None:
       print(" Error: fwidth is NOT given.")
    #
    shiftdist= (fwidth*j + i) * psize
    #
    fid.seek(shiftdist,0)
    outvalue = struct.unpack(fmt,fid.read(psize))
    #
    return outvalue
#    
def brij(infile,i,j,fwidth=None,psize=4,fmt='f'):
    #
    if fwidth is None:
       rsc = infile+'.rsc'
       info,ext = pSAR.roipac.rsc_read(rsc)
       fwidth = int(info['WIDTH'])
    #
    shiftdist= (fwidth*j + i) * psize
    #
    fid = open(infile,'rb')
    #
    # Find location of data
    #
    fid.seek(shiftdist,0)
    outvalue = struct.unpack(fmt,fid.read(psize))
    fid.close()
    #
    return outvalue
#
###############################################################################
def brll(infile,lon,lat,fwidth=None,psize=4,fmt='f'):
    #
    rsc  = infile+'.rsc'
    i,j  = pSAR.roipac.lonlat_to_index(rsc,[lon,lat])    
    return brij(infile,i,j,fwidth=None,psize=psize,fmt=fmt)
###############################################################################
def brlllist(infile,lllist,psize=4,fmt='f'):
    # 
    # updated by Wanpeng Feng, @Ottawa, 2016-12-27
    # lllist can be a numpy matrix
    #
    if isinstance(lllist,str):
       ll = np.loadtxt(lllist)
    else:
       ll = lllist
    #
    outdata  = np.zeros([ll.shape[0],3])
    ###########################################################################
    #
    rsc      = infile+'.rsc'
    info,ext = pSAR.roipac.rsc_read(rsc)
    i,j      = pSAR.roipac.ll_to_index_frominfo(info,ext,ll)
    #
    width    = int(float(info["WIDTH"]))
    length   = int(float(info["FILE_LENGTH"]))
    #
    fid      = open(infile,'rb')
    #
    i = i.astype('int')
    j = j.astype('int')
    #
    for index in range(ll.shape[0]):
        #
        outdata[index,0:2] = ll[index,0:2]
        #
        # i,j  = pSAR.roipac.lonlat_to_index(rsc,ll[index,:])
        #
        if (i[index] < width and j[index] < length):
           #
           outdata[index,2] = brij_id(fid,i[index],j[index],width,psize=psize,fmt=fmt)[0]
        else:
           outdata[index,2] = 9.
    #
    fid.close()
    return outdata
