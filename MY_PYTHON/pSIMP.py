#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 07:10:44 2017
It was decided to make an independent module for SIMP format that is defined in
PSOKINV saving fault geometry parameters. This format will be frequently used in
geodetic modelling. So it is good to have an independent python module. All functions 
in pDATA related to simp have been transferred to here.

@author: wanpeng Feng, @Ottawa, 2017-01-03
#
"""
###############################################################################   
import sys,numpy as np
import pSAR,math,cmath
try:
  import pGRD
  import pDATA
  import pGMT
except:
   print()
#
import os
import time

#
###############################################################################
def simp_pscmpinp(insimp,outinp,pts=np.zeros([180,2]),outpath='./output/',\
                  greenfunc=None,depth=0,lamd1=3.0e10,lamd2=3.0e10,\
                  days=[0]):
  '''
    create an inp configure file from a SIMP fault 
    outinp should be in format supported by PSGRN/PSCMP
    written by Wanpeng Feng, @SYSU, Guangzhou, 2020/04/05
  '''
  inp_header=\
'''#=======================================================================================================
# This is input file of FORTRAN program "pscmp2019" for modeling post-seismic deformation induced by
# earthquakes in multi-layered viscoelastic-gravitational media using the Green's function approach. The
# earthquke source is represented by an arbitrary number of rectangular dislocation planes.
#
# written by Rongjiang Wang
# Helmholtz Centre Potsdam
# GFZ German Research Centre for Geosciences
# e-mail: wang@gfz-potsdam.de
#
# Last modified: Potsdam, March, 2021
#
#################################################################
##                                                             ##
## Green's functions should have been prepared with the        ##
## program "psgrn2020" before the program "pscmp2020" is       ##
## started.                                                    ##
##                                                             ##
## For local Cartesian coordinate system, the Aki's convention ##
## is used, that is, x is northward, y is eastward, and z is   ##
## downward.                                                   ##
##                                                             ##
## If not specified otherwise, SI Unit System is used overall! ##
##                                                             ##
#################################################################
#=======================================================================================================
# OBSERVATION ARRAY
# =================
# 1. selection (0/1/2) for irregular observation positions (= 0) or a 1D observation profile (= 1) or a
#    rectangular 2D observation array (= 2): iposrec
#
#    IF (iposrec = 0 for irregular observation positions) THEN
#
# 2. number of positions: nrec
#
# 3. coordinates of the observations: (lat(i),lon(i)), i=1,nrec
#
#    ELSE IF (iposrec = 1 for regular 1D observation array) THEN
#
# 2. number of position samples of the profile: nrec
#
# 3. the start and end positions: (lat1,lon1), (lat2,lon2)
#
#    ELSE IF (iposrec = 2 for rectanglular 2D observation array) THEN
#
# 2. number of x samples, start and end values: nxrec, xrec1, xrec2
#
# 3. number of y samples, start and end values: nyrec, yrec1, yrec2
#
#    sequence of the positions in output data: lat(1),lon(1); ...; lat(nx),lon(1);
#    lat(1),lon(2); ...; lat(nx),lon(2); ...; lat(1),lon(ny); ...; lat(nx),lon(ny).
#=======================================================================================================
0
'''
  output_descript=\
'''#=======================================================================================================
# OUTPUTS
# =======
#
# 1. select (1/0) output for los displacement (only for snapshots, see below), x(n), y(e), and
#    z(d)-cosines of the unit vector from ground to satellite orbit: insar, xlos, ylos, zlos
#
#    if this option is selected, the snapshots will include additional data: Disp_LOS = los displacement
#    to the given satellite orbit.
#
# 2. select (1/0) output for Coulomb failure stress changes (only for snapshots, see below): icfs,
#    friction, Skempton ratio, strike, dip, and rake angles [deg] describing the uniform regional master
#    fault mechanism, the uniform regional principal stresses: sigma1, sigma2 and sigma3 [Pa] in
#    arbitrary order (the orietation of the pre-stress field will be derived by assuming that the master
#    fault is optimally oriented according to Coulomb failure criterion)
#
#    if this option is selected (icfs = 1), the snapshots will include additional data:
#
#    CFS_Max = increase of maximum CFS (can be oriented differently before and after the earthquake)
#    CFS_Mas = increase of CFS for the given master fault mechanism
#    CFS_Mas_Opt = increase of CFS on the master fault at the postseismic optimal rake direction
#    Sigma_Mas = increase of normal stress on the master fault
#    Rake_Mas_Opt = postseismic optimal rake on the master fault
#    CFS_Opt = increase of CFS for the postseismic optimal mechanism (in 3D space)
#    Sigma_Opt_1/2 = increase of normal stress on the first/second postseismic optimal faults
#    Strike_Opt_1/2, Dip_Opt_1/2, Rake_Opt_1/2 = the first/second postseismic optimal focal mechanisms
#
#    Note: the first 3D optimally orieted fault is the closest to the master fault.
#
# 3. output directory in char format: outdir
#
# 4. select outputs for displacement components (1/0 = yes/no): itout(i), i=1,3
#
# 5. the file names in char format for the x(n), y(e), and z(d) components: toutfile(i), i=1,3
#
# 6. select outputs for stress components (1/0 = yes/no): itout(i), i=4,9
#
# 7. the file names in char format for the xx(nn), yy(ee), zz(dd), xy(ne), yz(ed), and zx(dn)
#    components: toutfile(i), i=4,9
#
# 8. select outputs for vertical NS and EW tilt components, block rotation, geoid and gravity
#    changes (1/0 = yes/no): itout(i), i=10,14
#
# 9. the file names in char format for the NS tilt (positive if borehole top tilts to north), EW
#    tilt (positive if borehole top tilts to east), block rotation (clockwise positive), geoid and
#    gravity changes: toutfile(i), i=10,14
#
#    Note that all above outputs are time series with the time window as same as used for the Green's
#    functions
#
#10. number of scenario outputs ("snapshots": spatial distribution of all above observables at given
#    time points): nsc
#
#11. the time [day], and file name (in char format) for the 1. snapshot;
#12. the time [day], and file name (in char format) for the 2. snapshot;
#13. ...
#...
#=======================================================================================================
 0   -0.0068  0.3848  -0.9205      !cosines of LOS direction from ground to Evisat satellite orbit
 0     0.700  0.500  229.000   33.000  120.000    0.1E+07   -0.5E+07   -1.0E+07
'%s/'
  0                0               0
  'U_north.dat'    'U_east.dat'    'U_down.dat'
  0            0            0            0             0            0
  'S_nn.dat'   'S_ee.dat'   'S_dd.dat'   'S_ne.dat'    'S_ed.dat'   'S_dn.dat'
  0               0               0                0              0
  'Tilt_n.dat'    'Tilt_e.dat'    'Rotation.dat'   'geoid.dat'    'Gravity.dat'
  %d
'''
  greenfunc_path_descript=\
'''#=======================================================================================================
#
# GREEN'S FUNCTION DATABASE
# =========================
# 1. selection (0/1) of earth structure model (0 = homogeneous elastic halfspace, i.e., only for
#    co-seismic elastic Okada solutions, 1 = multi-layered viscoelastic halfspace, for which numerical
#    co- and post-seismic Green's functions calculated with psgrn2019 are needed): iesmodel
#
#    IF (iesmodel = 0 for analytical Okada solutions) THEN
#
# 2. observation depth [km] and the two Lame coefficients lambda and mue [Pa]
#
#    ELSE IF (iesmodel = 1 for numerical psgrn2019 Green's functions) THEN
#
# 2. directory where the Green's functions are stored: grndir
# 3. file names (without extensions!) for the 13 Green's functions:
#    3 displacement komponents (uz, ur, ut): green(i), i=1,3
#    6 stress components (szz, srr, stt, szr, srt, stz): green(i), i=4,9
#    radial and tangential components measured by a borehole tiltmeter,
#    rigid rotation around z-axis, geoid and gravity changes (tr, tt, rot, gd, gr):
#    green(i), i=10,14
#
#    Note that the extensions of the file names will be automatically considered.
#    They are ".ep", ".ss", ".ds" and ".cl" denoting the explosion (inflation)
#    strike-slip, the dip-slip and the compensated linear vector dipole sources,
#    respectively.
#
#=======================================================================================================
 %d
'%s'
%s
#======================================================================================================='''
  greenfunc_sub1=''' 'uz'  'ur'  'ut'
 'szz' 'srr' 'stt' 'szr' 'srt' 'stz'
 'tr'  'tt'  'rot' 'gd'  'gr' '''
  greenfunc_sub2='# No need to give green function files....'
  #
  fault_descript=\
'''#=======================================================================================================
# RECTANGULAR SUBFAULTS
# =====================
# 1. number of subfaults, latitude [deg] and east longitude [deg] of the regional reference point as
#    origin of the Cartesian coordinate system: ns, lat0, lon0
#
# 2. parameters for the 1. rectangular subfault: geographic coordinates (O_lat, O_lon) [deg] and O_depth
#    [km] of the local reference point on the present fault plane, length (along strike) [km] and width
#    (along down dip) [km], strike [deg], dip [deg], number of equi-size fault patches along the strike
#    (np_st) and along the dip (np_di) (total number of fault patches = np_st x np_di), and the start
#    time of the rupture; the following data lines describe the slip distribution on the present
#    sub-fault:
#
#    pos_s[km]  pos_d[km]  slip_stk[m]  slip_ddip[m]  open[m]
#
#    where (pos_s,pos_d) defines the position of the center of each patch in the local coordinate system
#    with the origin at the reference point: pos_s = distance along the length (positive in the strike
#    direction) pos_d = distance along the width (positive in the down-dip direction)
#    slip_stk = slip component in the strike direction
#    slip:ddip = slip component in the down dip direction (i.e., normal is positive and thrust is
#    negative)
#    open = fault opening (negative for closing)
#
#
# 3. ... for the 2. subfault ...
# ...
#                   N
#                  /
#                 /| strike
#                +------------------------
#                |\        p .            \ W
#                :-\      i .              \ i
#                |  \    l .                \ d
#                :90 \  S .                  \ t
#                |-dip\  .                    \ h
#                :     \. | rake               \
#                Z      -------------------------
#                              L e n g t h
#
#    Note that a point inflation can be simulated by three point openning faults (each causes a third
#    part of the volume of the point inflation) with orientation orthogonal to each other. the results
#    obtained should be multiplied by a scaling factor 3(1-nu)/(1+nu), where nu is the Poisson ratio at
#    the source. The scaling factor is the ratio of the seismic moment (energy) of an inflation source
#    to that of a tensile source inducing a plate openning with the same volume change.
#=======================================================================================================
# n_faults (Slip model by Ji Chen, USGS)
#-------------------------------------------------------------------------------
%d
'''
  faultpatch_descript=\
'''#-------------------------------------------------------------------------------------------------------
# n   O_lat   O_lon    O_depth length  width strike dip   np_st np_di start_time
# [-] [deg]   [deg]    [km]    [km]     [km] [deg]  [deg] [-]   [-]   [day]
#     pos_s   pos_d    slp_stk slp_ddip open
#     [km]    [km]     [m]     [m]      [m]
#-------------------------------------------------------------------------------------------------------
'''
  #
  fid = open(outinp,'w')
  fid.write(inp_header)
 
  # write points out to inp
  fid.write('%d\n' % pts.shape[0])
  #
  for i in range(pts.shape[0]):
    if (i+1) % 5 == 0:
        return_key = '\n'
    else:
        return_key = ''
    fid.write('(%f,%f) %s' % (pts[i,1],pts[i,0],return_key))
  #
  fid.write(output_descript % (outpath,len(days)))
  for cday in days:
      #
      fid.write(" %15.2f 'days_%d.dat'   |coseismic + postseismic of the first %d days\n" % (cday,cday,cday))
  #
  if greenfunc is None or not os.path.exists(greenfunc):
      flag_lay_func = 0
      greenfunc_str = '%5.3e %5.3e' % (lamd1,lamd2)
      substr = greenfunc_sub2
  else:
      flag_lay_func = 1
      greenfunc_str = greenfunc
      substr = greenfunc_sub1
  
  fid.write(greenfunc_path_descript % (flag_lay_func,greenfunc_str,substr))
  #
  # import simp in the coordiante of lonlat
  #
  simp,zone_num,zone_letter  = simp2pscmp(insimp)
  # simp,zone_num,zone_letter = import_simp(insimp)
  #
  #
  fid.write(fault_descript % simp.shape[0])
  fid.write(faultpatch_descript)
  
  for i in range(simp.shape[0]):
      fid.write('%-4d %f %f %f %f %f %f %f 1 1 0.00\n' % (i+1,simp[i,1],simp[i,0],\
                           simp[i,4], simp[i,6],simp[i,5],simp[i,2],simp[i,3]))
      fid.write('      %f %f %f %f %f\n' % (simp[i,6]/2,simp[i,5]/2,simp[i,7],simp[i,8]*-1,simp[i,9]))
  #
  fid.write('#==========================================end of input=================================================')
  fid.close()
  return True
#
def simp2pscmp(insimp):
    fpara,un,ul = import_simp(insimp)
    # convert from top-center to left-top
    tfpara = np.copy(fpara)
    for i in range(tfpara.shape[0]):
        #
        x0,y0,un,ul = pSAR.utm_conversion.from_latlon(tfpara[i,1],tfpara[i,0],\
                       force_zone_number=un,force_zone_letter=ul)
        #
        if tfpara[i,3]>90:
            tfpara[i,3] = 180 - tfpara[i,3]
            tfpara[i,2] = 180 + tfpara[i,2]
            tfpara[i,8] = tfpara[i,8] * -1
        #
        cfpara = np.copy(tfpara[i,:])
        cfpara[0] = x0/1000
        cfpara[1] = y0/1000
        cors,deps = simp_corners(cfpara) 
        lat,lon = pSAR.utm_conversion.to_latlon(cors[0,0]*1000,cors[0,1]*1000,\
                       un,zone_letter=ul)
        tfpara[i,0:2] = [lon,lat]
    return tfpara,un,ul
    #
def simp_enu2roi(inphs,ef,nf,uf,update=True):
    #
    ef_grd = ef+'.grd'
    nf_grd = nf+'.grd'
    uf_grd = uf+'.grd'
    #
    if not os.path.exists(ef_grd):
        pGMT.gmt_roi2grd(ef,ef_grd)
    if not os.path.exists(nf_grd):
        pGMT.gmt_roi2grd(nf,nf_grd)
    if not os.path.exists(uf_grd):
        pGMT.gmt_roi2grd(uf,uf_grd)
    #
    inphs_e_grd = inphs+'.e.grd'
    inphs_n_grd = inphs+'.n.grd'
    inphs_u_grd = inphs+'.u.grd'
    if not os.path.exists(inphs_e_grd) or update:
      pGMT.gmt_dem2roibyrsc(ef_grd,inphs+'.rsc',inphs_e_grd,gmt_r='-r')
    if not os.path.exists(inphs_n_grd) or update:
      pGMT.gmt_dem2roibyrsc(nf_grd,inphs+'.rsc',inphs_n_grd,gmt_r='-r')
    if not os.path.exists(inphs_u_grd) or update:
      pGMT.gmt_dem2roibyrsc(uf_grd,inphs+'.rsc',inphs_u_grd,gmt_r='-r')
    #
    e = pGRD.import_grd(inphs_e_grd)
    n = pGRD.import_grd(inphs_n_grd)
    u = pGRD.import_grd(inphs_u_grd)
    return e,n,u
    #
def sim_fparaconv(fpara,slocal,elocal):
    # 
    # this is a copy from psokinv (m codes)
    #
    # slocal and elocal can be
    #     0 -> up-center
    #     1 -> left-top
    #     2 -> rigth-top
    #     3 -> cc in the fault plane
    #     4 -> left-bottom
    #     5 -> right-bottom
    #    99 -> right-bottom, but with a depth on the top line
    #
    # elocal, only 0 is available, 2021/05/21
    # sim_fparaconv(fpara,99,0) should work for visco2pt5d
    #
    ofpara = np.copy(fpara)
    i = np.sqrt(-1+0j)
    #
    for ni in range(fpara.shape[0]):
        #
        strk = fpara[ni,2]
        dip  = fpara[ni,3]*np.pi/180.
        dep  = fpara[ni,4]
        wid  = fpara[ni,5]
        rlen = fpara[ni,6]
        x0   = fpara[ni,0]
        y0   = fpara[ni,1]
        #
        strkr = (90-strk) * np.pi/180.
        #
        if slocal == 0:
            ll = -0.5*rlen - i*wid*np.cos(dip)    # % low-left ,ul
            lr =  0.5*rlen - i*wid*np.cos(dip)    # % low-right,ur
            ul = -0.5*rlen + i*  0*np.cos(dip)    # % up-left  ,ll
            ur =  0.5*rlen + i*  0*np.cos(dip)    # % up-right ,lr
            sdep = wid*np.sin(dip)+dep
        elif slocal == 1:
            ll = 0    - i*wid*np.cos(dip)   #% low-left ,ul
            lr = rlen - i*wid*np.cos(dip)   #% low-right,ur
            ul = 0    + i*  0*np.cos(dip)   # % up-left  ,ll
            ur = rlen + i*  0*np.cos(dip)   #% up-right ,lr
            sdep= wid*np.sin(dip)+dep
        elif slocal == 2:
            ll = -rlen- i*wid*np.cos(dip)   #% low-left ,ul
            lr = 0    - i*wid*np.cos(dip)   #% low-right,ur
            ul = -rlen+ i*  0*np.cos(dip)   #% up-left  ,ll
            ur = 0    + i*  0*np.cos(dip)   #% up-right ,lr
            sdep= wid*np.sin(dip)+dep
        elif slocal == 99:
            #
            # source defined in visco2pt5d 
            #
            ll = -rlen  - i*0/2*np.cos(dip)    #% low-left ,ul
            lr =  0     - i*0/2*np.cos(dip)    #% low-right,ur
            ul =  -rlen + i*wid*np.cos(dip)    #% up-left  ,ll
            ur =  0     + i*wid*np.cos(dip)    #% up-right ,lr
            sdep = wid*np.sin(dip)+dep
            #
        elif slocal == 5: 
            # right-bottom
            #
            ll = -rlen  - i*0/2*np.cos(dip)    #% low-left ,ul
            lr =  0     - i*0/2*np.cos(dip)    #% low-right,ur
            ul =  -rlen + i*wid*np.cos(dip)    #% up-left  ,ll
            ur =  0     + i*wid*np.cos(dip)    #% up-right ,lr
            sdep = dep
        #
        ul = (x0+y0*i) + ul*np.exp(i*strkr)
        ur = (x0+y0*i) + ur*np.exp(i*strkr)
        ll = (x0+y0*i) + ll*np.exp(i*strkr)
        lr = (x0+y0*i) + lr*np.exp(i*strkr)
        #
        if elocal == 0:
            x = (np.real(ur)+np.real(ul))/2;
            y = (np.imag(ur)+np.imag(ul))/2;
            sdep = sdep-wid*np.sin(dip)
        # 
        ofpara[ni,0] = x
        ofpara[ni,1] = y
        ofpara[ni,4] = sdep
    return ofpara
#
###############################################################
def sim_fpara2rakes(fpara, slipthresh = 1.0):
    #
    rakes =  sim_fpara2rake(fpara)
    slip = np.sqrt(fpara[:,7]**2 + fpara[:,8]**2)
    return rakes,[np.max(rakes[slip>=slipthresh]),\
                  np.min(rakes[slip>=slipthresh])]
#    
def sim_fpara2rake(fpara):
    #
    if len(fpara) == 10:
        fpara = np.reshape(fpara,[1,10])
    return np.arctan2(fpara[:,8],fpara[:,7])*180./np.pi
#
###############################################################################
#
def sim_mfpara2gmt(fpara,outgmt_root,mode='xh',rot=None):
    #    
    #
    polygon_gmt = outgmt_root+'.poly.gmt'
    psvelo_gmt  = outgmt_root+'.psvelo.gmt'
    #
    fid_polygon = open(polygon_gmt,'w')
    fid_psvelo  = open(psvelo_gmt,'w')
    #
    for i in range(fpara.shape[0]):
        #
        outdata,tolslip,vel = simp_fpara4gmt(fpara[i,:],mode=mode,rot=rot)
        #
        fid_polygon.write('>-Z%f\n' % tolslip)
        #
        for k in range(outdata.shape[0]):
            fid_polygon.write('%f %f \n' % (outdata[k,0],outdata[k,1]))
        #
        fid_psvelo.write('%f %f %f %f %f %f %f\n' % (vel[0],vel[1],vel[2],\
                                                     vel[3],vel[4],vel[5],\
                                                     vel[6]))
    #
    fid_polygon.close()
    fid_psvelo.close()
    return polygon_gmt,psvelo_gmt
#
###############################################################################
#
def simp_fpara4gmt(fpara,mode='xh',rot=None):
    '''
    rot is a flag showing if fpara should be rotated based on the strike angle
        None is default for not rotating, otherwise [refx,refy,strike] 
        should be given.
    mode flag for output coordinates, xh for x and depth
         xy, for x and y
         yh, for y and depth
    '''
    fpara = np.reshape(fpara,10)
    xy,depths = simp_corners(fpara)
    if rot is not None:
       rx,ry = simp_rotax(xy[:,0],xy[:,1],rot[2],rot[0],rot[1])
    else:
       rx,ry = np.copy(xy[:,0]),np.copy(xy[:,1])
    #
    if mode.upper() == "XY":
        outdata = np.vstack((rx,ry)).T
    if mode.upper() == "XH":
        outdata = np.vstack((rx,depths)).T
    if mode.upper() == 'YH':
        outdata = np.vstack((ry,depths)).T
    #
    tolslip = np.sqrt(fpara[7]**2+fpara[8]**2)
    # print(' FWP: ***** %d' % outdata.shape[0])
    vel     = [np.mean(outdata[0:4,0]),np.mean(outdata[0:4,1]),fpara[7],fpara[8]*-1.,\
                       0,0,0]
    return outdata,tolslip,vel
#   
###############################################################################
#     
def simp_momentalongstrike(insimp,refll=None):
    '''
       xory flag for selecting the dimension to separate the slip model
       by Wanpeng Feng, @NRcan, 2017-12-07
    '''
    #
    faults,zn,zl = import_simp(insimp)
    if refll is None:
        refx,refy = 0.,0.
    else:
        refx,refy,a,b = pSAR.utm_conversion.from_latlon(refll[1],refll[0],\
                                        force_zone_number=zn,\
                                        force_zone_letter=zl)
        refx,refy = refx/1000.,refy/1000.
    #
    rotx,roty = simp_rotax(faults[:,0],faults[:,1],np.mean(faults[:,2]),\
                           refx,refy)
    #
    urotx = np.unique(rotx)
    urotx = np.sort(urotx)
    momentprof = np.zeros([urotx.shape[0],2])
    #
    # mstep = pSAR.util.meanstep(urotx)
    udepths = np.unique(faults[:,4])
    cfpara = faults[faults[:,4]==udepths[0],:]
    mstep = (np.max(urotx)-np.min(urotx))/cfpara.shape[0]
    #
    nfault = faults.shape[0]/cfpara.shape[0]
    print(" pSIMP: %d faults at the samle depth" % cfpara.shape[0])
    print(" pSIMP: %d faults at the same x-post" % nfault)
    for i in range(urotx.shape[0]):
        #
        flag1 = rotx>(urotx[i]-mstep/5.)
        flag2 = rotx<(urotx[i]+mstep/5.)
        #
        flag = flag1 & flag2
        cfparas = faults[flag,:]
        # print("Strike VS moment with %d faults" % cfparas.shape[0])
        a,b = simp_moment(cfparas)
        momentprof[i,:] = [urotx[i],b]
    return momentprof
#
###############################################################################
#    
def simp_momentprofile(insimp):
    #
    faults,zn,zl = import_simp(insimp)
    depths = faults[:,4]
    udepths = np.unique(depths)
    udepths = np.sort(udepths)
    momentprof = np.zeros([udepths.shape[0],2])
    for i in range(udepths.shape[0]):
        #
        cfparas = faults[faults[:,4]==udepths[i],:]
        a,b = simp_moment(cfparas)
        momentprof[i,:] = [b,udepths[i]]
    return momentprof
#
def export_simp(fpara,outsimp,zone=None,nl=None):
    simp_fault2simp(fpara,outsimp,zone=zone,nl=nl);
    return True
#
def simp_fault2simp(fault,outsimp,zone=None,nl=None):
    #
    if len(fault.shape) == 1:
        fault = np.reshape(fault,[1,10])
    #
    if zone is None or nl is None:
        print(" pSIMP: zone and nl have been provided...")
        return False
    #
    nl = nl.upper()
    with open(outsimp,'w') as fid:
      fid.write("# Updated on %s\n" % time.strftime("%Y%m%dT%H:%M:%S"))
      fid.write("# UTM ZONE: %s%s\n" % (str(zone),str(nl)))
      fid.write("# Number of faults: %d \n" % fault.shape[0])
      fid.write("# x(km) y(km) str(deg) dip(deg) dep(km) wid(km) leng(km) s_slip(m) d_slip(m) o_slip(m)\n")
      for i in range(fault.shape[0]):
          fid.write('%f %f %f %f %f %f %f %f %f %f\n' % \
            (fault[i,0],fault[i,1],fault[i,2],fault[i,3],fault[i,4],\
             fault[i,5],fault[i,6],fault[i,7],fault[i,8],fault[i,9]))
    #
    return True
    #
def simp_azi(p1,p2):
    #
    p1,p2 = np.array(p1),np.array(p2)
    x,y = p2[0] - p1[0], p2[1] - p1[1]
    #
    ang1 = np.arctan2(x,y)
    #
    return np.rad2deg( ang1 % (2 * np.pi))

def simp_moment2mag(mom):
    return np.log10(mom) * 2./3. - 6.033    
#
def simp_mag2moment(mag):
    return 10. ** ((mag+6.033)*(3./2.))
#    
def simp4mntratealongstrike(simpfile,npoints=20,refx=None,refy=None,mu=3.2e10):
    #
    x,y = simp4momentrate(simpfile,refx=refx,refy=refy,mu=mu)
    xpos= (np.max(x)-np.min(x))/npoints * 1.5
    #
    tmpx = np.linspace(np.min(x)-xpos,np.max(x)+xpos,npoints+2)
    tmpy = np.copy(tmpx) * 0
    #
    for ck in range(tmpx.shape[0]):
      index = ck+1
      if index <= tmpx.shape[0]-1:
        x0 = tmpx[index] 
        x1 = tmpx[index] + xpos
        flag1 = x >= x0
        flag2 = x <  x1
        flag  = np.logical_and(flag1,flag2)
        tmpy[index] = np.sum(y[flag])
        #
    return tmpx,tmpy
#
###############################################################################
def simp2moment(simp,mu=3.2e10):
    fpara,un,l = import_simp(simp)
    sm,tolm = simp_moment(fpara,mu=3.2e10)
    mw = simp_moment2mag(tolm)
    return tolm, mw
#
###############################################################################
#
def simp_moment(simp,mu=3.2e10):
    #
    smoment   = mu * simp[:,5] * simp[:,6] * \
                (np.sqrt(simp[:,7]**2 + simp[:,8]**2)) * 10**6
    return smoment,np.sum(smoment)
#
#
def simp4momentrate(simpfile,refy=None,refx=None,mu=3.2e10):
#
    simp,zone_num,zone_letter = import_simp(simpfile)
    refkx,refky,n,l           = pSAR.utm_conversion.from_latlon(refy,refx,\
                force_zone_number=zone_num,force_zone_letter=zone_letter)
    #
    #
    refkx   = refkx / 1000
    refky   = refky / 1000
    simp_x  = np.copy(simp[:,0])
    simp_y  = np.copy(simp[:,1])
    simp_d  = np.copy(simp_x) * 0
    #
    for ni in range(simp.shape[0]):
       #
       polys,depth = simp_corners(simp[ni,:])
       simp_x[ni]  = np.mean(polys[0:4,0])
       simp_y[ni]  = np.mean(polys[0:4,1])
       simp_d[ni]  = np.mean(depth[0:4])
    #
    newx,newy = simp_rotax(simp_x,simp_y,simp[0,2],refkx,refky)
    smoment   = mu * simp[:,5] * simp[:,6] * (np.sqrt(simp[:,7]**2 + \
                simp[:,8]**2)) * 10**6
    return newx,smoment
#    
###############################################################################   
#
def simp4vec(simpfile,refx=None,refy=None,ax2='d',ax3='syn'):
    #
    simp,outdata= simp4xyzc(simpfile,None,isout=False,refx=refx,\
                            refy=refy,ax2=ax2,ax3=ax3)
    slipvec     = np.zeros([simp.shape[0],7])
    slipvec[:,0]= outdata[:,0]
    slipvec[:,1]= outdata[:,1]
    slipvec[:,2]= simp[:,7]
    slipvec[:,3]= simp[:,8]
    slipvec[:,6]= 0.1
    return slipvec
#
def simp4xyzc(simpfile,outname,isout=True,refx=None,refy=None,ax2='d',ax3='syn'):
    #
    if refx is None:
       #
       print("")
       print("simp4xyzc(simpfile,outname,refx=None,refy=None,ax2='d',ax3='syn')")
       print(" return xyz data from a simp faul model and xy will be the center of patches")
       print(" refx and refy are necessary for now...")
       sys.exit(-1)
    #
    simp,zone_num,zone_letter = import_simp(simpfile)
    refkx,refky,n,l           = pSAR.utm_conversion.from_latlon(refy,refx,\
                                force_zone_number=zone_num,\
                                force_zone_letter=zone_letter)
    #
    #
    refkx   = refkx / 1000
    refky   = refky / 1000
    simp_x  = np.copy(simp[:,0])
    simp_y  = np.copy(simp[:,1])
    simp_d  = np.copy(simp_x) * 0
    #
    for ni in range(simp.shape[0]):
       #
       polys,depth = simp_corners(simp[ni,:])
       simp_x[ni]  = np.mean(polys[0:4,0])
       simp_y[ni]  = np.mean(polys[0:4,1])
       simp_d[ni]  = np.mean(depth[0:4])
       #print(depth.shape)
    #
    ###################################################
    #
    newx,newy = simp_rotax(simp_x,simp_y,simp[0,2],refkx,refky)
    if ax2.upper()=="D":
       pts = np.vstack((newx,simp_d)).T
    if ax2.upper()=="Y":
        pts = np.vstack((newx,simp_y)).T
    #
    if ax3.upper()=="SYN":
        ax3input = np.sqrt(simp[:,7]**2 + simp[:,8]**2)
    if ax3.upper()=="STR":
        ax3input = simp[:,7]
    if ax3.upper()=="DIP":
        ax3input = simp[:,8]
    if ax3.upper()=="OPEN":
        ax3input = simp[:,9]
    #
    outdata = np.hstack((pts,np.reshape(ax3input,[ax3input.shape[0],1])))
    if isout:
       np.savetxt(outname,outdata,fmt='%f %f %f',newline='\n')
    #
    return simp,outdata
#        
###############################################################################
#
def import_simp(simpfmodel):
    #
    # read distributed fault model in simp format
    # we can control fault numbers by reading the num directly from the file..
    # updated by FWP, @SYSU, Zhuhai, 2021/01/26
    #
    fid = open(simpfmodel,'r')
    counter = 0
    mfault  = []
    gofault = False
    fnum = None
    fcounter = 0
    #
    for cline in fid:
        counter += 1
        if '# Number of faults' in cline:
           tmp = cline.split('\n')[0].split(':')
           fnum = int(tmp[1].split()[0])
        if '# UTM ZONE' in cline:
           zone = cline.split('\n')[0]
           zone = zone.split(':')
           zone = ''.join(zone[1:])
           zoneNum = int(zone[0:-1])
           zoneL   = zone[-1]
        #
        if gofault:
           tmp = cline.split('\n')[0]
           tmp = tmp.split()
           #
           # updated by FWP, @SYSU, Zhuhai, 2021/01/16
           #
           if fcounter < fnum:
             mfault.append(pDATA.strarr_to_num(tmp)) 
             fcounter += 1
        #
        if 'dip(deg) dep(km) wid(km)' in cline: 
           gofault = True
    #
    fid.close()
    mfault = np.array(mfault) 
    # print(zoneNum)
    return mfault,zoneNum,zoneL 
###############################################################################
def simp2geoxyz(insimp):
    #
    utmxyz,zonenum,zoneletter = simp2xyz(insimp)
    #
    outdata = np.copy(utmxyz)
    #
    for ni in range(utmxyz.shape[0]):
        outdata[ni,1],outdata[ni,0] = pSAR.utm_conversion.to_latlon(\
               utmxyz[ni,0]*1000, utmxyz[ni,1]*1000,\
               zonenum,zone_letter=zoneletter)
    #                  
    return outdata,zonenum,zoneletter
    #
###############################################################################   
def simp2grd(insimp,outgrd,refll=None,model='xz',psize=None,gmt_r='',sf=1):
    #
    xyz,uz,ul = simp2xyz(insimp)
    deps      = xyz[:,2]
    udeps     = np.unique(deps)
    ny        = udeps.shape[0]
    nx        = xyz[xyz[:,2]==udeps[0],:].shape[0]
    #
    if refll is None:
       # 
       refx,refy = 0.,0.
    else:
       refx,refy,a,b = pSAR.utm_conversion.from_latlon(refll[1],refll[0],\
                                                       force_zone_number=uz,\
                                                       force_zone_letter=ul)
       refx,refy = refx/1000.,refy/1000.
    #
    ux,uy = simp_rotax(xyz[:,0],xyz[:,1],np.mean(xyz[:,3]),refx,refy)
    npt   = ux.shape[0]
    #
    slip = np.sqrt(xyz[:,4]**2 + xyz[:,5]**2)
    #
    if model.upper() == 'XZ':
       #
       outdata = np.hstack((np.reshape(ux,[npt,1]),\
                                       np.reshape(xyz[:,2],[npt,1]),\
                                       np.reshape(slip,[npt,1])))
       minx,maxx = np.min(ux),np.max(ux)
       miny,maxy = np.min(xyz[:,2]),np.max(xyz[:,2])
       #
    if model.upper() == 'XY':
       #
       outdata = np.hstack((np.reshape(ux,[npt,1]),\
                            np.reshape(uy,[npt,1]),\
                            np.reshape(slip,[npt,1])))
       #
       minx,maxx = np.min(ux),np.max(ux)
       miny,maxy = np.min(uy),np.max(uy)
    #
    gmt_ext = '-R%f/%f/%f/%f' % (minx,maxx,miny,maxy)
    # pGMT
    if psize is None:
        xsize = (maxx-minx)/(nx-1)
        ysize = (maxy-miny)/(ny-1)
    else:
        xsize = psize[0]
        ysize = psize[1]
    #
    xsize = xsize / sf
    ysize = ysize / sf
    #
    gmt_inc = '-I%f/%f' % (xsize,ysize)
    outxyz = os.path.basename(insimp)+'.xyz'
    np.savetxt(outxyz,outdata,fmt='%f %f %f',newline='\n')
    #
    flag = pGMT.gmt_xyz2grd(outxyz,outgrd,gext=gmt_ext,surface=True,\
                gmt_r=gmt_r,ginc=gmt_inc,binary=False,A=1,di='-di0')
    #
    return flag,xsize,ysize
#
###############################################################################
#
def simp2xyz(insimp):
    #
    fpara,zonenum,zoneLetter = import_simp(insimp)
    outdata = []
    #
    for ni in range(fpara.shape[0]):
        utmpoly,utmdepth = simp_corners(fpara[ni,:])
        cx = utmpoly[0:4,0]
        cy = utmpoly[0:4,1]
        cz = utmdepth[0:4]
        outdata.append([np.mean(cx),np.mean(cy),\
                        np.mean(cz),fpara[ni,2],fpara[ni,7],fpara[ni,8],fpara[ni,9]])
    #
    outdata = np.array(outdata)                    
    return outdata,zonenum,zoneLetter
#
###############################################################################  
#
def simp_rotaxll(x,y,strike,cx0,cy0):
    #
    x0, y0, zn, zl = pSAR.utm_conversion.from_latlon(cy0,cx0)
    x0, y0 = x0 / 1000., y0 / 1000.
    x, y = pSAR.utm_conversion.lltoutm(y,x,force_zone_number=zn,force_zone_letter=zl)
    x, y = x / 1000., y / 1000.
    outx, outy = simp_rotax(x,y,strike,x0,y0)
    return outx, outy
#
###############################################################################   
#def simp_topline(fpara):
#    # can only be applied to a single fault
#    #
#    if fpara[4] == 0:
#       utmpoly,utmdepth = simp_corners(fpara)
#       return utmpoly,utmdepth
#    else:
#       return None
def simp_rotax(x,y,strike,cx0,cy0):
    #
    #
    x, y    = np.array(x), np.array(y)
    PI      = np.pi
    strkr   = (90.-strike)*PI/180.;
    istrkr  = complex(0,strkr)
    
    xystart = complex(cx0,cy0);
    cx      = np.copy(x)
    cy      = np.copy(y)
    #
    # print(cx)
    #
    for ni in range(cx.shape[0]):
       outdata = (complex(x[ni],y[ni])-xystart) * cmath.exp(-1. * istrkr)
       cx[ni] = outdata.real
       cy[ni] = outdata.imag   
    return cx,cy   
###############################################################################  
def simp_to_maskcorner(insimp,refax=1):
    #
    # refax 1 for y , 0 for x
    #
    # utmxyz,zonenum,zoneletter = simp2xyz(insimp)   
    outdata,zonenum,zoneletter = simp2geoxyz(insimp)
    zdepth = outdata[:,2]
    cpolytop  = outdata[zdepth==np.min(zdepth),:]
    cpolybot  = outdata[zdepth==np.max(zdepth),:]
    #
    p1 = cpolytop[cpolytop[:,refax] == np.min(cpolytop[:,refax]),:]
    p2 = cpolytop[cpolytop[:,refax] == np.max(cpolytop[:,refax]),:]
    p3 = cpolybot[cpolybot[:,refax] == np.max(cpolybot[:,refax]),:] 
    p4 = cpolybot[cpolybot[:,refax] == np.min(cpolybot[:,refax]),:]
    #
    outp = np.zeros([5,3])
    #
    # print(p1.shape)
    #
    outp[0,:] = p1[0,0:3]
    outp[1,:] = p2[0,0:3]
    outp[2,:] = p3[0,0:3]
    outp[3,:] = p4[0,0:3]
    outp[4,:] = p1[0,0:3]
    #
    return outp
###############################################################################
def simp2geoext(insimp):
    #
    outll,outdeps = simp2geopoly(insimp)
    #
    lats = []
    lons = []
    for i in range(outll.shape[0]):
        #
        cpoly = np.reshape(outll[i,:],[2,5])
        lats.append(cpoly[1,:])
        lons.append(cpoly[0,:])
    #
    lats = np.array(lats)
    lons = np.array(lons)
    return [np.min(lons.ravel()),np.max(lons.ravel()),\
            np.min(lats.ravel()),np.max(lats.ravel())],\
            [np.min(outdeps.ravel()),np.max(outdeps.ravel())]
    #
def simp_fpara2geopoly(fpara,zn,zl):
    #
    fpara = np.array(fpara)
    if len(fpara.shape) == 1 and len(fpara)==10:
        fpara = np.reshape(fpara,(1,10))
    #
    outll = []
    outdepth = []
    #
    for ni in range(fpara.shape[0]):
        #
        # cslip = math.sqrt(fpara[ni,7]**2 + fpara[ni,8]**2)
        #
        utmpoly,utmdepth = simp_corners(fpara[ni,:])
        lon,lat = pSAR.utm_conversion.utmtoll(utmpoly[:,0] * 1000,\
                                              utmpoly[:,1] * 1000,\
                                              zn,\
                                              force_zone_letter=zl)
        outll.append(np.hstack([lon,lat]))
        outdepth.append(utmdepth.ravel())
    #
    outll = np.array(outll)
    outdepth = np.array(outdepth)
    #
    return outll,outdepth
    #
def simp2geopoly(insimp):
    #
    #
    fpara,zonenum,zoneLetter = import_simp(insimp)
    # 
    outll = []
    outdepth = []
    #
    for ni in range(fpara.shape[0]):
        #
        # cslip = math.sqrt(fpara[ni,7]**2 + fpara[ni,8]**2)
        #
        utmpoly,utmdepth = simp_corners(fpara[ni,:])
        lon,lat = pSAR.utm_conversion.utmtoll(utmpoly[:,0] * 1000,\
                                              utmpoly[:,1] * 1000,\
                                              zonenum,\
                                              force_zone_letter=zoneLetter)
        outll.append(np.hstack([lon,lat]))
        outdepth.append(utmdepth.ravel())
    #
    outll = np.array(outll)
    outdepth = np.array(outdepth)
    #
    return outll,outdepth

###############################################################################
def simp2gmtpoly(insimp,outgmt,isopenning=False,dep=False):
    #
    #
    fpara,zonenum,zoneLetter = import_simp(insimp)
    #
    outid = open(outgmt,'w')
    for ni in range(fpara.shape[0]):
        #
        if not isopenning:
           cslip = math.sqrt(fpara[ni,7]**2 + fpara[ni,8]**2)
        else:
           cslip = fpara[ni,9]
        utmpoly,utmdepth = simp_corners(fpara[ni,:])
        #
        lon,lat = pDATA.utm2lonlat(utmpoly[:,0] * 1000,utmpoly[:,1] * 1000,\
                                   zonenum,zoneLetter)
        outid.write("> -Z%f \n" % cslip)
        for nk in range(5):
          if dep:
            outid.write("%f %f %f\n" % (lon[nk],lat[nk],utmdepth[nk]))
          else:
            outid.write("%f %f\n" % (lon[nk],lat[nk]))
    #
    outid.close()
    return None
###############################################################################
def simp_2pts2cutm(lons,lats,zone=None,nl=None):
    if zone is None or nl is None:
       x,y,zone_num,zone_letter = pSAR.utm_conversion.from_latlon(\
                                            np.mean(lats),np.mean(lons))
    else:
       x,y,zone_num,zone_letter = pSAR.utm_conversion.from_latlon(\
                                            np.mean(lats),np.mean(lons),\
                                           force_zone_number=int(zone),\
                                           force_zone_letter=nl) 
    return x,y,zone_num,zone_letter
    #
def simp_2pts2length(lons,lats,zone=None,nl=None):
    '''
    '''
    x0,y0,un,ul = simp_2pts2cutm(lons,lats,zone=None,nl=None)
    x,y = pSAR.utm_conversion.lltoutm(lats,lons,\
                                force_zone_number=un,force_zone_letter=ul)
    x,y = x/1000.,y/1000.
    return np.sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2)
    #
#
def simp_2pts2fault(lons,lats,dip=88.,wid=10,dep=0.,zone=None,nl=None):
    '''
    zone utm zone number
    nl   utm zone letter
    dip  dip angle of fault
    dep  burial depth
    '''
    x0,y0,zone,nl = simp_2pts2cutm(lons,lats,zone=zone,nl=nl)
    strike        = simp_2ptstostrike(lons,lats,zone=zone,nl=nl)
    if strike < 0:
        strike = 360.0 + strike
    #
    flen          = simp_2pts2length(lons,lats,zone=zone,nl=nl)
    outfault      = np.zeros([1,10])
    outfault[0,0] = x0/1000.
    outfault[0,1] = y0/1000.
    outfault[0,2] = strike
    outfault[0,3] = dip
    outfault[0,4] = dep
    outfault[0,5] = wid
    outfault[0,6] = flen
    return outfault,zone,nl
    #
def simp_2ptstostrike(lons,lats,zone=None,nl=None):
    #
    x0,y0,zone_num,zone_letter = simp_2pts2cutm(lons,lats,zone=zone,nl=nl)
    x1,y1,z,l = pSAR.utm_conversion.from_latlon(lats[0],lons[0],\
                                                force_zone_number=zone_num,\
                                                force_zone_letter=zone_letter)
    x2,y2,z,l = pSAR.utm_conversion.from_latlon(lats[1],lons[1],\
                                                force_zone_number=zone_num,\
                                                force_zone_letter=zone_letter)
    #
    dx = (x2-x1) / 1000.
    dy = (y2-y1) / 1000.
    return 90 - np.arctan2(dy,dx) * 180 / np.pi 
#
###############################################################################
#
def simp_pt2profile(lon,lat,strike,rlen,mode="l"):
    #
    #
    x,y,zone_num,zone_letter = pSAR.utm_conversion.from_latlon(lat,lon)
    #
    simp = np.array([x/1000,y/1000,strike,89.9999,0,10,rlen*2,0,0,0])
    #
    polys,rdepth = simp_corners(simp)
    #
    if mode.upper()=="L":
       x1,y1        = polys[1,0],polys[1,1]
       outvalue     = pSAR.utm_conversion.to_latlon(x1*1000.,y1*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       outprofile   = np.array([[lon,lat],[outvalue[1],outvalue[0]]])
    if mode.upper()=="R":
       x1,y1        = polys[0,0],polys[0,1]
       outvalue     = pSAR.utm_conversion.to_latlon(x1*1000.,y1*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       outprofile   = np.array([[lon,lat],[outvalue[1],outvalue[0]]])
    if mode.upper()=="C":
       x1,y1        = polys[0,0],polys[0,1]
       outvalue1    = pSAR.utm_conversion.to_latlon(x1*1000.,y1*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       x1,y1        = polys[1,0],polys[1,1]
       outvalue2    = pSAR.utm_conversion.to_latlon(x1*1000.,y1*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       outprofile   = np.array([[outvalue1[1],outvalue1[0]],\
                                [outvalue2[1],outvalue2[0]]])
    if mode.upper()=="EXT":
       x1,y1        = polys[1,0],polys[1,1]
       outvalue1    = pSAR.utm_conversion.to_latlon(x1*1000.,y1*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       x2,y2        = polys[0,0],polys[0,1]
       outvalue2    = pSAR.utm_conversion.to_latlon(x2*1000.,y2*1000,\
                        zone_number=zone_num,zone_letter=zone_letter)
       outprofile   = np.array([[outvalue1[1],outvalue1[0]],\
                                [outvalue2[1],outvalue2[0]]])
    #
    return outprofile

###############################################################################
#
def refll_to_profile(ll,strike,plen=50,spacing=10,isutm=False):
    #
    if not isutm:
      ux,uy,zone_num,zone_letter = pSAR.utm_conversion.from_latlon(ll[1],ll[0])
    else:
      ux,uy = ll[0]*1000.,ll[1]*1000.
      zone_num,zone_letter = None,None
    #
    simp      = np.array([ux/1000,uy/1000,strike,90,0,10,plen,0,0,0])
    poly,deps = simp_corners(simp)
    ###########################################################################
    #
    # UR
    #
    x0   = poly[0,0]
    y0   = poly[0,1]
    # UL
    x1   = poly[1,0]
    y1   = poly[1,1]
    flen = np.sqrt((x0-x1)**2 + (y0-y1)**2)
    # print(x0,y0,x1,y1,flen)
    fnum = int(flen / spacing)
    kmx  = np.linspace(x0,x1,num=fnum) * 1000
    kmy  = np.linspace(y0,y1,num=fnum) * 1000
    lons,lats = pDATA.utm2lonlat(kmx,kmy,zone_num,zone_letter)
    return lons,lats,zone_num,zone_letter
###############################################################################
#
def import_simp_inlonlat(simpfmodel):
    # read simp fault model to numpy array in lonlat
    #
    mfault,zoneNum,zoneLetter = import_simp(simpfmodel)
    mfault = simp_to_lonlat(mfault,zoneNum,zoneLetter)
    return mfault
###############################################################################
#
def simp_roifromlist(inlist):
    #
    #
    datalist = []
    counter  = 0
    with open(inlist,'r') as fid:
        for cline in fid:
            cline = pSAR.util.bytestoutf8(cline)
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
###############################################################################
def simp_fpara2profile(fpara,step):
    #
    polys,depths = simp_corners(fpara)
    x0,x1 = polys[0,0],polys[1,0]
    y0,y1 = polys[0,1],polys[1,1]
    npts = int(fpara[6] / step)
    if npts < 1:
       npts = 1
    mx,my = np.linspace(x0,x1,num = npts+1),np.linspace(y0,y1,num=npts+1)
    return mx,my
#    
###############################################################################
def simp_to_lonlat(simpfault,zoneNum,zoneLetter):
    #
    for ni in range(simpfault.shape[0]):
        x = simpfault[ni,0] * 1000
        y = simpfault[ni,1] * 1000
        lat,lon = pSAR.utm_conversion.to_latlon(x,y,zoneNum,zoneLetter)
        simpfault[ni,0] = lon
        simpfault[ni,1] = lat
    #
    return simpfault

###############################################################################
def simp_topline(simp,u=None,l=None):
    #
    fpara = simp_fpara2whole(simp,0)
    poly,deps = simp_corners(fpara)
    poly,deps = poly[0:3,:],deps[0:3]
    cpoly = poly[deps == 0,:]
    if u is not None:
       lpoly = np.copy(cpoly)
       lats,lons = pSAR.utm_conversion.utmtoll(cpoly[:,0]*1000,cpoly[:,1]*1000,\
                                             u,force_zone_letter=l)
       lpoly[:,0],lpoly[:,1] = lons,lats
       cpoly = lpoly
    return cpoly
#
###############################################################################
def simp_updatefpara(simp,strike=None,dip=None,model='foc'):
    #
    csimp = np.copy(simp)
    if len(csimp)>10:
        print(" Error: only single band of fault model is required")
        sys.exit(-1)
    #
    csimp = np.reshape(csimp,10)
    xy,deps = simp_corners(simp)
    #
    if model.upper() == 'TOP':
        newloc = [csimp[0],csimp[1],csimp[4]]
    if model.upper() == 'FOC':
       newloc = simp_fpara2center(simp)
    # 
    csimp[0] = newloc[0]
    csimp[1] = newloc[1]
    csimp[4] = newloc[2]
    #
    if strike is not None:
        csimp[2] = strike
    if dip is not None:
        csimp[3] = dip
    #
    nsimp = simp_fpara2whole(csimp,np.min(deps),maxdepth=np.max(deps))
    nsimp = np.reshape(nsimp,[1,10])
    #
    return nsimp
#
def simp_mfpara2center(mfpara):
    '''
    To return the centers of fault models
    
    '''
    outdata = np.zeros([mfpara.shape[0],3])
    for ni in range(mfpara.shape[0]):
        cdata = simp_fpara2center(mfpara[ni,:])
        outdata[ni,:] = cdata
    #
    return outdata
#
def simp_fpara2center(simp):
#
    #
    xy,deps = simp_corners(simp)
    meanx   = np.mean(xy[0:4,0])
    meany   = np.mean(xy[0:4,1])
    meand   = np.mean(deps[0:4])
    #
    return [meanx,meany,meand]
#
def simp_fpara2whole(simp,rdepth,maxdepth=None):
    '''
    Calculate a new fault plane based on the given <simp>. The new plane should 
    be in the same plane, but constrained by a depth range, [rdepth,maxdepth]
    
    by Wanpeng Feng, @Ottawa, 2017-11-23
    '''
    simp  = np.reshape(simp,10)
    rsimp = np.copy(simp)
    rsimp = rsimp
    #
    DD     = np.pi/180.
    dip    = simp[3] * DD
    xcent  = simp[0]
    ycent  = simp[1]
    rlen   = simp[6]
    cwidth = simp[5]
    strike = simp[2]
    depth  = simp[4]
    #
    if maxdepth == None:
       maxdepth = np.sin(dip)*cwidth+ depth  
    #
    if maxdepth < rdepth:
       print(" ERROR: maxdepth(%f) is smaller than rdepth(%f)" % \
             (maxdepth,rdepth))
    #
    ul     = complex(-0.5*rlen, 0)  
    strkr  = (90-strike) * DD 
    istrkr = complex(0,strkr)
    ul     = complex(xcent, ycent) + ul*cmath.exp(istrkr)
    top    = depth
    #
    rwidth = top/np.sin(dip)
    rul    = complex(0*rlen,   rwidth*np.cos(dip))   #  % up right corner
    rur    = complex(  rlen,   rwidth*np.cos(dip))   #  % up right corner
    rul    = ul + rul*np.exp(istrkr)
    rur    = ul + rur*np.exp(istrkr)
    #
    rsimp[0]  = (rul.real + rur.real)/2.
    rsimp[1]  = (rul.imag + rur.imag)/2.
    rsimp[5]  = (maxdepth - rdepth)/np.sin(dip)
    #
    top = rdepth
    #
    topc        = complex(0,   -1. * top/np.sin(dip) * np.cos(dip))
    topc        = complex(rsimp[0], rsimp[1]) + topc * cmath.exp(istrkr)
    #
    rsimp[0]    = topc.real
    rsimp[1]    = topc.imag
    rsimp[4]    = rdepth
    return rsimp
#
def simp_llcorners(simp,u,l):
    #
    poly,depths = simp_corners(simp)
    llpoly = np.copy(poly) * 0.
    poly   = poly * 1000.
    lats,lons = pSAR.utm_conversion.utmtoll(poly[:,0],\
                                                poly[:,1],\
                                                u,force_zone_letter=l)
    llpoly[:,0] = lons
    llpoly[:,1] = lats
    return llpoly
#
def simp_fpara2whole_batch(insimp,depth,maxdepth=None):
    #
    outsimp = np.copy(insimp)
    #
    for i in range(outsimp.shape[0]):
        outsimp[i,:] = simp_fpara2whole(insimp[i,:],depth,maxdepth=maxdepth)
    #
    return outsimp
#
def simp_corners(simp):
    '''
    order of ouput: ul,ur,lr,ll,ul 
    '''
    # Note that this can be only applied to a single fault
    #
    simp  = np.reshape(simp,[10])
    polys = []
    zdepth= []
    # calculate corners of a rectangle fault plane
    PI     = 3.14159265
    rlen   = simp[6]
    wid    = simp[5]
    dip    = simp[3] / 180. * PI
    dep    = simp[4]
    strike = simp[2]
    x0     = simp[0]
    y0     = simp[1]
    ###########################################################################
    ll = complex(-0.5*rlen, -1*wid*math.cos(dip))  # low-left ,ul
    lr = complex( 0.5*rlen, -1*wid*math.cos(dip))  # low-right,ur
    ul = complex(-0.5*rlen,      0*math.cos(dip))  # up-left  ,ll
    ur = complex( 0.5*rlen,      0*math.cos(dip))  # up-right ,lr
    strkr = (90.-strike)*PI/180.;
    istrkr= complex(0,strkr)
    ul = complex(x0,y0) + ul*cmath.exp(istrkr)
    ur = complex(x0,y0) + ur*cmath.exp(istrkr)
    ll = complex(x0,y0) + ll*cmath.exp(istrkr)
    lr = complex(x0,y0) + lr*cmath.exp(istrkr)
    # UL
    polys.append([ul.real,ul.imag])  # UL
    polys.append([ur.real,ur.imag])  # UR
    polys.append([lr.real,lr.imag])  # LR
    polys.append([ll.real,ll.imag])  # LL
    polys.append([ul.real,ul.imag])  # UL for a colosure
    polys = np.array(polys)
    # Depth
    zdepth.append(dep)  # UL
    zdepth.append(dep)  # UR
    zdepth.append(dep+wid*math.sin(dip)) # LR
    zdepth.append(dep+wid*math.sin(dip)) # LR
    zdepth.append(dep)                   # UL
    zdepth = np.array(zdepth)
    #
    return polys,zdepth
