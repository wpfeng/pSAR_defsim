#!/usr/bin/env python
#
#
'''
Created on Fri Nov 11 17:28:24 2022

create time series
#
@author: wafeng
'''
#
import pokada
import pSAR
import pSIMP
import numpy as np
import datetime
import os
import sys
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import glob
from scipy.fftpack import fftshift, fft2, ifft2
###
#
def rednoise(N1, N2, r):
    if N1 % 2:
        M1 = N1 + 1
    else:
        M1 = N1
    #
    if N2 % 2:
        M2 = N2 + 1
    else:
        M2 = N2
    #
    x = np.random.randn(M1, M2)
    X = fftshift(fft2(x))
    n1 = np.concatenate((np.arange(-M1//2, 0), np.arange(1, M1//2 + 1)))
    n2 = np.concatenate((np.arange(-M2//2, 0), np.arange(1, M2//2 + 1)))
    n1, n2 = np.meshgrid(n1, n2)
    kk = np.sqrt(n1**2 + n2**2)**r
    X = X / kk
    y = ifft2(fftshift(X))
    y = np.real(y[:N1, :N2])
    #
    y = y - np.mean(y)
    yrms = np.sqrt(np.mean(y**2))
    y = y / yrms
    #
    #y = y * std_deviation / 3 
    return y
    #
#
def create_los(los_dir, bins,emag,vx,vy,vz):
    #
    if not os.path.exists(los_dir):
        #
        os.makedirs(los_dir)
        #
    #
    infs = glob.glob(inf_dir+'/geo*.phs')
    if len(infs)>0:
        print(' %d interferograms have been processed... ' % len(infs))
        #
    else:
      for i in range(len(bins)-1):
        #
        #
        bin_1_e = bins[i]
        bin_1_n = bins[i].replace('_E.bin','_N.bin')
        bin_1_u = bins[i].replace('_E.bin','_U.bin')
        #
        stime = os.path.basename(bin_1_e).split('_E.bin')[0]
        outlos = '%s/%s.los' % (los_dir,stime)
        #
        if not os.path.exists(outlos):
          de = pSAR.roipac.roi_read(bin_1_e)
          dn = pSAR.roipac.roi_read(bin_1_n)
          du = pSAR.roipac.roi_read(bin_1_u)
          #
          info,ext = pSAR.roipac.rsc_read(bins[i]+'.rsc')
          #
          aps_err = rednoise(de.shape[0],de.shape[1],r) * emag / 3
          dlos =  de * vx + dn*vy + du*vz + aps_err    #(np.random.rand(de.shape[0],de.shape[1])-0.5)* emag * 2
          #
          #
          #
          pSAR.roipac.roi_write(dlos,outlos)
          #
          pSAR.roipac.info_to_rsc(info,outlos+'.rsc')

def create_infs(inf_dir,bins):
    #
    if not os.path.exists(inf_dir):
        #
        os.makedirs(inf_dir)
        #
    #
    infs = glob.glob(inf_dir+'/geo*.phs')
    if len(infs)>0:
        print(' %d interferograms have been processed... ' % len(infs))
        #
    else:
      for i in range(len(bins)-1):
        #
        info,ext = pSAR.roipac.rsc_read(bins[i]+'.rsc')
        #
        bin_1_e = bins[i]
        #
        dlos1 = pSAR.roipac.roi_read(bin_1_e)
        #
        # data1 =  (de * vx + dn*vy + du*vz + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag)*-1*np.pi/wavelength
        data1 =  dlos1 *-4*np.pi/wavelength
        #
        mdate = bin_1_e.split('/')[1].split('T')[0]
        #
        for j in range(len(bins)):
            #
            bin_2_e = bins[j]
            #
            dlos2 = pSAR.roipac.roi_read(bin_2_e)
            #
            # data2 =  (de * vx + dn*vy + du*vz + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag*2)*-1*np.pi/wavelength
            data2 =  dlos2*-4*np.pi/wavelength
            #
            #
            sdate = bin_2_e.split('/')[1].split('T')[0]
            #
            diffdays = pSAR.ts.diff2dates(mdate,sdate)
            season_flag = abs(diffdays - np.round(diffdays/365)*365)
            #
            #
            if diffdays>0:
               #
               if diffdays < 96 or season_flag < 30:

                 # data2 = pSAR.roipac.roi_read(bin_2_e)
                 info2,_ = pSAR.roipac.rsc_read(bin_2_e+'.rsc')
                 info['STIME'] = info2['MTIME']
                 info['SLAVE'] = info2['MASTER']
                 info['PERPB']  = float(info2['PERPB']) - float(info['PERPB'])
                 #
                 outdata = data2 - data1
                 #
                 #
                 outinf = '%s/geo_%s_%s_%s.phs' % (inf_dir,mdate,sdate,trackinfo)
                 #
                 pSAR.roipac.roi_write(outdata,outinf)
                 #
                 # info['PERB'] = 0
                 pSAR.roipac.info_to_rsc(info,outinf+'.rsc')

def losvec(inc,azi,mode='rng',ldir='right'):
    #
    inc_rad = inc * np.pi / 180
    azi_rad = azi * np.pi / 180
    if ldir.upper()=="LEFT":
      azi_rad = azi * np.pi / 180 + 3.14159265
    #
    n_rng =      np.sin(azi_rad) * np.sin(inc_rad)
    e_rng = -1 * np.cos(azi_rad) * np.sin(inc_rad)
    u_rng =      np.cos(inc_rad)
    #
    if mode.upper() == 'RNG':
        return e_rng,n_rng,u_rng
    else:
      #
      n_azi = math.cos(azi_rad)
      e_azi = math.sin(azi_rad)
      u_azi = 0
      return e_azi,n_azi,u_azi
#
def dis_xy_locs(fpara,pixelsize,ext=None):
    #
    polys,_ = pSIMP.simp_corners(fpara)
    cext = [polys[0:4,0].min(),polys[0:4,0].max(),polys[0:4,1].min(),polys[0:4,1].max()]
    #
    if ext is None:
       xwid = cext[1]-cext[0]
       ywid = cext[3]-cext[2]
       ext = [cext[0]-xwid*3.75,cext[1]+xwid*3.75, \
              cext[2]-ywid*3.75,cext[3]+ywid*3.75]
    #
    intx = int((ext[1]-ext[0])/pixelsize)
    inty = int((ext[3]-ext[2])/pixelsize)
    #
    x = np.linspace(ext[0],ext[1],intx)
    y = np.linspace(ext[3],ext[2],inty)
    #
    ux,uy = np.meshgrid(x,y)
    #
    return ux,uy
#
def slipshistory(jds,cosd_jd,log_para,cosslip,linearate):
    #
    islips = (jds-jds[0]) * linearate/365.
    #
    slips = np.copy(islips)*0
    if cosd_jd > jds[-1]:
        slips[:] = 0
        #
    elif cosd_jd < jds[0]:
        tjds = jds - cosd_jd
        slips[:] =  log_para[0]*np.log(1+tjds/log_para[1])
    else:
        #
        tjds = jds - cosd_jd
        #
        slips[tjds>=0] = cosslip + log_para[0]*np.log(1+tjds[tjds>0]/log_para[1])
        slips[tjds<0]  = 0
    #
    return slips, islips
#
def syn2insar(outsyn,yyyymmdd,samples):
    #
    outinsar = []
    # 
    for i in range(yyyymmdd.shape[0]-1):
        #
        for j in range(1,samples+1):
            #
            # print('%d-%d' % (yyyymmdd[i],yyyymmdd[j+i]))
            #
            if (j+i) < yyyymmdd.shape[0]:
               # print('%d-%d' % (yyyymmdd[i],yyyymmdd[j+i]))
               clos = outsyn[j+i] - outsyn[i]
               outinsar.append([yyyymmdd[i],yyyymmdd[j+i],clos])
    #
    return np.array(outinsar)
#
def toDT(date0,date1,timestep):
    #
    d0 = datetime.datetime.strptime(str(date0),'%Y%m%d')
    dts = []
    yyyymmdd = []
    #
    flag = 1
    while flag>0:
        #
        #
        cdt = d0 + datetime.timedelta(days=timestep)
        #
        cdate = cdt.strftime('%Y%m%d')
        #
        if int(cdate)<=date1:
            dts.append(cdt)
            yyyymmdd.append(cdate)
            d0 = cdt
        else:
            flag = 0
        
        #
    #
    return np.array(dts),np.array(yyyymmdd).astype('int')
#
def dis_ts(intime,v0):
    #
    return intime*v0/365
#
if len(sys.argv)<3:
    #
    helpstr = \
        '''
        %s <start_time> <end_time> -hhmmss [01:01:01] -chhmmss [01:01:01] -cosd [20990101] 
                                   -infdir [InSAR] 
                                   -disdir [Forward_DIS_DIR]
                                   -timestep [12 in default]
                                   -fpara [0,0,80,70,2,5,5,0,0,0 in default]
                                   -log_para [3.5, 50]
                                   -linearate [0.005]
                                   -cosslip   [0.6]
                                   -ext [None]
                                   -wavelength [0.0556 in default]
                                   -epi [-117.5,33.5]
                                   -inc [35 in default]
                                   -azi [-12 in default]
                                   -ldir [right or left, right in default]
                                   -mode [rng or azi, rng in default]
                                   -emag [0.001m in default, error magnitude]
                                   -trackinfo [T001 in default]
                                   -pixelsize [1 in default]

        Since 2023/04/14, the script was developed to simulate surface deformation due to 
        variable faulting behaviors in the elastic half space Earth model.

        The code has been released along with the study below:
        
        Zhang, Z., Feng, W., Xu, X., & Samsonov, S. (2023). Performance of 
        Common Scene Stacking Atmospheric Correction on Nonlinear InSAR 
        Deformation Retrieval. Remote Sensing, 15(22), 5399.

        '''
    #
    print(helpstr % sys.argv[0])
    #
    sys.exit(-1)
    #
#
pSAR.util.log(sys.argv)
logfile = os.path.basename(sys.argv[0]).split('.py')[0]+'.log'
#
#
intime          =  int(sys.argv[1])
endtime         =  int(sys.argv[2])
hhmmss          = '01:01:01'
chhmmss         = '01:01:01'
cosd            =  20990101
inf_dir         = 'InSAR'
out_forward_dir = 'Forward_DIS_DIR'
timestep        = 12
log_para        = [0.05,50]
linearate       = 0.005
cosslip         = 0.85
fpara           = [0,0,80,70,2,5,5,0,0,0]
ext             = None
pixelsize       = 1
epi             = [-117.5,33.5]
wavelength      = 0.0556
plot            = False
inc             = 35
azi             = -12
ldir            = 'right'
emag            = 0.0001
trackinfo       = 'T001'
r               = 1.8         # The recommended value range is 4/3 - 8/3.
# 
# InSAR viewing mode, rng or azi
mode            = 'rng'
#
#
#
#
#
for i,key in enumerate(sys.argv):
    #
    if key == '-trackinfo':
       trackinfo = sys.argv[i+1]
    if key == '-emag':
        emag = float(sys.argv[i+1]) 
    if key == '-azi':
       azi = float(sys.argv[i+1])
    if key == '-inc':
       inc = float(sys.argv[i+1])
    if key == '-epi':
       epi = [float(cloc) for cloc in sys.argv[i+1].split(',')]
    if key == '-wavelength':
       wavelength = float(sys.argv[i+1])
    if key == '-plot':
       plot = True
    if key == '-ext':
       ext = [float(cext) for cext in sys.argv[i+1].split(',')]
    if key == '-pixelsize':
       pixelsize = float(sys.argv[i+1])
    if key == '-cosslip':
       cosslip = float(sys.argv[i+1])
    if key == '-linearate':
       linearate = float(sys.argv[i+1])
       #
    if key == '-log_para':
       log_para = [float(cpara) for cpara in sys.argv[i+1].split(',')]
    if key == '-fpara':
       fpara = [float(cloc) for cloc in sys.argv[i+1].split(',')]
    if key == '-timestep':
       timestep = int(sys.argv[i+1])
    if key == '-hhmmss':
        hhmmss = sys.argv[i+1]
    if key == '-chhmmss': 
        chhmmss = sys.argv[i+1]
    if key == '-cosd':
        cosd = int(sys.argv[i+1])
    if key == '-infdir':
        inf_dir = sys.argv[i+1]
    if key == '-disdir':
        out_forward_dir = sys.argv[i+1]
    if key == '-r':
        r   =  sys.argv[i+1]
# 
if True:
    #
    #
    if abs(azi)<30:
        dirstr = 'ASC'
    else:
        dirstr = 'DES'
    #
    fpara = np.array(fpara)
    #
    samples = 100 
    #
    # cosdate = '20160617'
    # dates = np.array([i for i in range(0,1000,12)])
    #
    dts,yyyymmdds = toDT(intime,endtime,timestep)
    #
    jds = np.array([pSAR.ts.yyyymmdd2jd(cdate,hhmmss=hhmmss,fmt='%Y%m%dT%H:%M:%S') for cdate in yyyymmdds])
    cosd_jd = pSAR.ts.yyyymmdd2jd(cosd,hhmmss=chhmmss)
    #
    #
    # slips = 1.8*np.log(1+jds/50)
    #
    slips, islips = slipshistory(jds,cosd_jd,log_para,cosslip,linearate)
    #
    if plot:
      plt.plot(jds-cosd_jd,slips,'-or')
      plt.show()
      #
    #
    out_forward_dir = 'Forward_DIS_DIR'
    #
    if not os.path.exists(out_forward_dir):
        #
        os.makedirs(out_forward_dir)
        #
    #
    mmerr = emag*1000
    err_str = 'MM_%03.1f' % mmerr 
    los_dir_err = 'LOS_DIR_ERR_%s' % err_str
    inf_dir = 'InSAR_%s_%s' % (err_str,dirstr)
    inf_dir_noerr = 'InSAR_%s' % dirstr
    #
    #
    ux,uy = dis_xy_locs(fpara,pixelsize,ext=ext)
    #
    print(" *** ux and uy are ready for use, with %d and %d pixels..." % (ux.shape[0],ux.shape[1]))
    #
    if plot:
      plt.subplot(121)
      plt.imshow(ux)
      plt.colorbar()
      plt.subplot(122)
      plt.imshow(uy)
      plt.colorbar()
      plt.show()
      #
      sys.exit(-1)
    #
    x_step = pixelsize
    y_step = pixelsize*-1
    ext = [ux.ravel().min(),ux.ravel().max(),uy.ravel().min(),uy.ravel().max()]
    #
    x,y,un,ul = pSAR.utm_conversion.from_latlon(epi[1],epi[0])
    x0 = x
    y0 = y
    #
    utmext = [ext[0]+x/1000,ext[1]+x/1000,ext[2]+y/1000,ext[3]+y/1000]
    #
    latmin,lonmin = pSAR.utm_conversion.utmtoll(np.array([utmext[0]*1000]),np.array([utmext[2]*1000]),un,force_zone_letter=ul)
    latmax,lonmax = pSAR.utm_conversion.utmtoll(np.array([utmext[1]*1000]),np.array([utmext[3]*1000]),un,force_zone_letter=ul)
    #
    latmin,lonmin = latmin[0],lonmin[0]
    latmax,lonmax = latmax[0],lonmax[0]
    #
    geoext = [lonmin,lonmax,latmin,latmax]
    x_step = (lonmax-lonmin)/(ux.shape[1]-1)
    y_step = (latmin-latmax)/(ux.shape[0]-1)
    #
    #
    alp = 0.25
    #
    info,_ = pSAR.roipac.roipac_info()
    #
    hhmmss_str = ''.join(hhmmss.split(':'))
    #
    #
    polys,depths = pSIMP.simp_corners(fpara)
    #
    # cfpara, an interseismic slipping fault
    # edited by Wanpeng Feng, 2023/11/26, to make sure interseismic model can be modified based on inputs...
    #
    cfpara = [polys[2:4,0].mean(),polys[2:4,1].mean(),fpara[2],fpara[3],depths[2],fpara[5],fpara[6],0,0,0]
    #
    fparas = np.zeros([2,10])
    fparas[0,:] = fpara
    fparas[1,:] = cfpara
    fparas[:,7] = 0.5
    #
    llfpara = np.copy(fparas)
    llfpara[:,0] = llfpara[:,0]+x0/1000
    llfpara[:,1] = llfpara[:,1]+y0/1000
    #
    #
    if not os.path.exists('fparas_XY.simp'):
      #
      pSIMP.simp_fault2simp(fparas,'fparas_XY.simp',zone=un,nl=ul)
      pSIMP.simp_fault2simp(llfpara,'fparas_UTM.simp',zone=un,nl=ul)
      os.system('pSAR_simp2llsimp.py fparas_UTM.simp fparas_LL.simp')
    #
    #
    if plot:
      plt.plot(slips)
      plt.plot(islips)
      plt.show()
    #
    ## 
    # calculate state vectors...
    #
    with open(logfile,'a') as fid:
        fid.write('###++++++++++++++++++++++++++++++++++++++++++++++\n')
        fid.write('# inc:%f, azi:%f, mode:%s\n' % (inc,azi,mode))
        fid.write('###++++++++++++++++++++++++++++++++++++++++++++++\n')
    #
    vx,vy,vz = losvec(inc,azi,ldir=ldir,mode=mode)
    print(vx,vy,vz)
    #
    for i,slip in enumerate(slips):
      #
      # formating output name...
      #
      outname = out_forward_dir+'/%dT%s.bin' % (yyyymmdds[i],hhmmss_str)
      outname_e = out_forward_dir+'/%dT%s_E.bin' % (yyyymmdds[i],hhmmss_str)
      outname_n = out_forward_dir+'/%dT%s_N.bin' % (yyyymmdds[i],hhmmss_str)
      outname_u = out_forward_dir+'/%dT%s_U.bin' % (yyyymmdds[i],hhmmss_str)
      #
      #
      if not os.path.exists(outname_e):
        #
        fparas[0,7] = slip 
        fparas[1,7] = islips[i]
        #
        dx,dy,dz = pokada.dc2dfromfaults(fparas,ux,uy,alpha=alp,disp=True)
        dx_err = dx + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag*2
        dy_err = dy + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag*2
        dz_err = dz + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag*2
        #
        # dlos = (dx * vx + dy*vy + dz*vz + (np.random.rand(ux.shape[0],ux.shape[1])-0.5)*emag)*-1*np.pi/wavelength
        #
        #
        if not os.path.exists(outname_e):
          print(' ... %s is ready...' % outname_e)
          pSAR.roipac.roi_write(dx,outname_e,outrsc=True)
          pSAR.roipac.roi_write(dy,outname_n,outrsc=True)
          pSAR.roipac.roi_write(dz,outname_u,outrsc=True)
        #
        #
        info['FILE_LENGTH'] = dz.shape[0]
        info['WIDTH']       = dz.shape[1]
        info['X_FIRST']     = geoext[0] 
        info['Y_FIRST']     = geoext[3]
        info['X_STEP']      = x_step 
        info['Y_STEP']      = y_step
        info['WAVELENGTH']  = wavelength
        info['HEADING_DEG'] = azi
        info['INCIDENCE']   = inc
        #
        yyyy = int(str(yyyymmdds[i])[0:4])
        mm   = int(str(yyyymmdds[i])[4:6])
        dd   = int(str(yyyymmdds[i])[6:8])
        #
        info['MTIME']       = '%d-%d-%dT%s.01Z' % (yyyy,mm,dd,hhmmss)
        info['MASTER']      = yyyymmdds[i]
        info['PERPB']       = (np.random.rand() - 0.5) * 400 
        #
        #
        #
        pSAR.roipac.info_to_rsc(info,outname_e+'.rsc')
        pSAR.roipac.info_to_rsc(info,outname_n+'.rsc')
        pSAR.roipac.info_to_rsc(info,outname_u+'.rsc')
        #
    #
    bins = glob.glob('%s/*_E.bin' % out_forward_dir)
    bins.sort()
    #
    # creat los + errors
    # which an be applied to create infs and extract error levels...
    #
    create_los(los_dir_err, bins,emag,vx,vy,vz)
    #
    los_dir_noerr = 'LOS_DIR'
    create_los(los_dir_noerr, bins,0,vx,vy,vz)
    #
    #
    #
    losfiles = glob.glob('%s/*.los' % los_dir_err)
    losfiles.sort()
    print(" ... now let us make interferograms")
    create_infs(inf_dir,losfiles)
    #
