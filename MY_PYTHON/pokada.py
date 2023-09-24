#!/usr/bin/env python
#

"""
  Created on Fri Sep  30 05:08:43 2016
  Updated on Sun June 06 06:06:00 2021
     Since this version, the point source for explosive source is available...
  #
  This is an integrated module of elastic dislocation inlcuding both 2d and 3d versions. 
  The two dislocation functions, okadaDC3D and okadaDC2D, are compultely based on the original  
  Fortran codes released with  classic half-space elatic dislocation journal papers:
  #
    Okada, Y., 1992. BSSA.
    Okada, Y., 1985. BSSA.
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                
  The kernel elsatic disloation codes are only a transmission from .f to .py. But this python version
  has been modified to adapt matrix computation in use, which can allow one to calculate 
  surface deformation significantly faster. This seems very important to search fault geometry parameters 
  in a very considerably loose parameter space. My previous Multipeak Particle Swarm 
  Optimization (MPSO) algorithm is also being prepared under python. A new version global geodetic
  inversion package will also be coming soon.
  #
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Note that the original "return code" RECT was currently out of use. Then here may be some Nan
  value returned when locations are just on the fault trace. So take care of your results carefully.
  In fact, original fortran functions have also sigularity issues. Some potential solutions can be
  found in the below webpage
  http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html
  
  This has included Okada2D and Okada3D codes. I just used my previous Matlab 
  based okada model to validate this module. Two results from the both versions
  are completely identical.
  If you find it useful and put it input your own package, please acknowlege our 
  efforts correctly by citing below papers. We should release this after they are 
  submitted. 
  <1>
  <2>
  Any error reports should always be appreciated. Please report them to Wanpeng Feng,
  wanpeng.feng@hotmail.com
  
  
@author: Wanpeng Feng
"""
import numpy as np
import numpy.matlib as npm
import numba
#
###############################################################################
# point source
def okadadc3d0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4):
#      SUBROUTINE  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,
#     *               UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
#      IMPLICIT REAL*8 (A-H,O-Z)
#      REAL*4   ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,
#     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
#C
#C********************************************************************
#C*****                                                          *****
#C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****
#C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
#C*****                         CODED BY  Y.OKADA ... SEP.1991   *****
#C*****                         REVISED     NOV.1991, MAY.2002   *****
#C*****                                                          *****
#C********************************************************************
#C
#C***** INPUT
#C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
#C*****   X,Y,Z : COORDINATE OF OBSERVING POINT
#C*****   DEPTH : SOURCE DEPTH
#C*****   DIP   : DIP-ANGLE (DEGREE)
#C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
#C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2
#C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3
#C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4
#C
#C***** OUTPUT
#C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /
#C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )
#C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /
#C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )
#C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE
#C*****   IRET        : RETURN CODE
#C*****               :   =0....NORMAL
#C*****               :   =1....SINGULAR
#C*****               :   =2....POSITIVE Z WAS GIVEN
#C
#      COMMON /C1/DUMMY(8),R
#      DIMENSION  U(12),DUA(12),DUB(12),DUC(12)

#      DATA  F0/0.D0/
      F0 = 0.0
#C-----
      IRET=0
      if (Z > 0.):
        Z = Z * -1
      #
#C-----
      #for I in range(12):
      U   = [np.copy(X)*F0] * 12
      DUA = [np.copy(X)*F0] * 12
      DUB = [np.copy(X)*F0] * 12
      DUC = [np.copy(X)*F0] * 12
      #  111 CONTINUE
      AALPHA = ALPHA
      DDIP   = DIP
      ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D = \
             DCCON0(AALPHA,DDIP)
#C======================================
#C=====  REAL-SOURCE CONTRIBUTION  =====
#C======================================
      XX=X
      YY=Y
      ZZ=Z
      DD=DEPTH+Z
      # calling DCCON1
      P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ = \
         DCCON1(XX,YY,DD,SD,CD)
      if (R == F0):
        return U
      #
      PP1=POT1
      PP2=POT2
      PP3=POT3
      PP4=POT4
      #
      DUA = UA0(XX,YY,DD,PP1,PP2,PP3,PP4,\
                ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ)
# C-----#
      for I in range(12):
        if (I < 9):
           U[I]=U[I]-DUA[I]
        if (I >= 9):
           U[I]=U[I]+DUA[I]
#
#C=======================================
#C=====  IMAGE-SOURCE CONTRIBUTION  =====
#C=======================================
      DD=DEPTH-Z
      # calling DCCON1
      P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ = \
            DCCON1(XX,YY,DD,SD,CD)
      DUA = UA0(XX,YY,DD,PP1,PP2,PP3,PP4,\
                ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ)
      DUB = UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,\
               ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D, \
             P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,\
             UY,VY,WY,UZ,VZ,WZ)
      DUC = UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,\
                ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D, \
             P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3)
#C-----
      for I in range(12):
        DU=DUA[I]+DUB[I]+ZZ*DUC[I]
        if (I>=9):
           DU=DU+DUC[I-9]
        U[I]=U[I]+DU
#
      return U
#      UX=U[0]
#      UY=U[1]
#     UZ=U[2]
#      UXX=U[3]
#      UYX=U[4]
#      UZX=U[5]
#      UXY=U[6]
#      UYY=U[7]
#      UZY=U[8]
#      UXZ=U[9]
#      UYZ=U[10]
#      UZZ=U[11]
#      
#      return UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
############################################################################
def UA0(X,Y,D,POT1,POT2,POT3,POT4,\
        ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
        P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ):
#      SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U)
#      IMPLICIT REAL*8 (A-H,O-Z)
#      DIMENSION U(12),DU(12)
#C
#C********************************************************************
#C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****
#C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
#C********************************************************************
#C
#C***** INPUT
#C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM
#C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
#C***** OUTPUT
#C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
#C
#      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
#      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,
#     *           UY,VY,WY,UZ,VZ,WZ
#      DATA F0,F1,F3/0.D0,1.D0,3.D0/
#      DATA PI2/6.283185307179586D0/
      #
      F0 = 0.0
      F1 = 1.0
      F3 = 3.0
      PI2 = np.pi*2     
      #
      U = [np.copy(X)*F0]*12
      DU = [np.copy(X)*F0] * 12
#C-----
#C======================================
#C=====  STRIKE-SLIP CONTRIBUTION  =====
#C======================================
      if (POT1 != F0):
        DU[0]  = ALP1*Q/R3    +ALP2*X2*QR
        DU[1]  = ALP1*X/R3*SD +ALP2*XY*QR
        DU[2]  =-ALP1*X/R3*CD +ALP2*X*D*QR
        DU[3]  = X*QR*(-ALP1 +ALP2*(F1+A5) )
        DU[4]  = ALP1*A3/R3*SD +ALP2*Y*QR*A5
        DU[5]  =-ALP1*A3/R3*CD +ALP2*D*QR*A5
        DU[6]  = ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY
        DU[7]  = F3*X/R5*(-ALP1*Y*SD +ALP2*(Y*UY+Q) )
        DU[8]  = F3*X/R5*( ALP1*Y*CD +ALP2*D*UY )
        DU[9]  = ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ
        DU[10] = F3*X/R5*( ALP1*D*SD +ALP2*Y*UZ )
        DU[11] = F3*X/R5*(-ALP1*D*CD +ALP2*(D*UZ-Q) )
        for I in range(12):
           U[I]=U[I]+POT1/PI2*DU[I]
#      ENDIF
#C===================================
#C=====  DIP-SLIP CONTRIBUTION  =====
#C===================================
      if (POT2 != F0):
        DU[0]  =            ALP2*X*P*QR
        DU[1]  = ALP1*S/R3 +ALP2*Y*P*QR
        DU[2]  =-ALP1*T/R3 +ALP2*D*P*QR
        DU[3]  =                 ALP2*P*QR*A5
        DU[4]  =-ALP1*F3*X*S/R5 -ALP2*Y*P*QRX
        DU[5]  = ALP1*F3*X*T/R5 -ALP2*D*P*QRX
        DU[6]  =                          ALP2*F3*X/R5*VY
        DU[7]  = ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR)
        DU[8]  =-ALP1*(C2D/R3-F3*Y*T/R5) +ALP2*F3*D/R5*VY
        DU[9]  =                          ALP2*F3*X/R5*VZ
        DU[10] = ALP1*(C2D/R3+F3*D*S/R5) +ALP2*F3*Y/R5*VZ
        DU[11]= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR)
        for I in range(12):
           U[I]=U[I]+POT2/PI2*DU[I]
#      ENDIF
#C========================================
#C=====  TENSILE-FAULT CONTRIBUTION  =====
#C========================================
      if (POT3 != F0):
        DU[0]  = ALP1*X/R3 -ALP2*X*Q*QR
        DU[1]  = ALP1*T/R3 -ALP2*Y*Q*QR
        DU[2]  = ALP1*S/R3 -ALP2*D*Q*QR
        DU[3]  = ALP1*A3/R3     -ALP2*Q*QR*A5
        DU[4]  =-ALP1*F3*X*T/R5 +ALP2*Y*Q*QRX
        DU[5]  =-ALP1*F3*X*S/R5 +ALP2*D*Q*QRX
        DU[6]  =-ALP1*F3*XY/R5           -ALP2*X*QR*WY
        DU[7]  = ALP1*(C2D/R3-F3*Y*T/R5) -ALP2*(Y*WY+Q)*QR
        DU[8]  = ALP1*(S2D/R3-F3*Y*S/R5) -ALP2*D*QR*WY
        DU[9]  = ALP1*F3*X*D/R5          -ALP2*X*QR*WZ
        DU[10] =-ALP1*(S2D/R3-F3*D*T/R5) -ALP2*Y*QR*WZ
        DU[11] = ALP1*(C2D/R3+F3*D*S/R5) -ALP2*(D*WZ-Q)*QR
        for I in range(12):
           U[I]=U[I]+POT3/PI2*DU[I]
#      ENDIF
#C=========================================
#C=====  INFLATE SOURCE CONTRIBUTION  =====
#C=========================================
      if (POT4 != F0):
        DU[0]  =-ALP1*X/R3
        DU[1]  =-ALP1*Y/R3
        DU[2]  =-ALP1*D/R3
        DU[3]  =-ALP1*A3/R3
        DU[4]  = ALP1*F3*XY/R5
        DU[5]  = ALP1*F3*X*D/R5
        DU[6]  = DU[4]
        DU[7]  =-ALP1*B3/R3
        DU[8]  = ALP1*F3*Y*D/R5
        DU[9]  =-DU[5]
        DU[10] =-DU[8]
        DU[11] = ALP1*C3/R3
        for I in range(12):
           U[I]=U[I]+POT4/PI2*DU[I]
      return U
#
def UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,\
        ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D, \
             P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,\
             UY,VY,WY,UZ,VZ,WZ):
#      SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)
#      IMPLICIT REAL*8 (A-H,O-Z)
#      DIMENSION U(12),DU(12)
#C
#C********************************************************************
#C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
#C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
#C********************************************************************
#C
#C***** INPUT
#C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM
#C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
#C***** OUTPUT
#C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
#C
#      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
#      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,
#     *           UY,VY,WY,UZ,VZ,WZ
#      DATA F0,F1,F2,F3,F4,F5,F8,F9
#     *        /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,8.D0,9.D0/
#      DATA PI2/6.283185307179586D0/
#C-----
      #
      F0 = 0
      F1 = 1.
      F2 = 2.0
      F3 = 3.0
      F4 = 4.0
      F5 = 5.0
      F8 = 8.0
      F9 = 9.0
      PI2 = np.pi*2
      #
      U = [np.copy(X)*F0] * 12
      DU = [np.copy(X)*F0] * 12
      C=D+Z
      RD=R+D
      D12=F1/(R*RD*RD)
      D32=D12*(F2*R+D)/R2
      D33=D12*(F3*R+D)/(R2*RD)
      D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)
      D54=D12*(F5*R2+F4*R*D+D2)/R3*D12
# C-----
      FI1= Y*(D12-X2*D33)
      FI2= X*(D12-Y2*D33)
      FI3= X/R3-FI2
      FI4=-XY*D32
      FI5= F1/(R*RD)-X2*D32
      FJ1=-F3*XY*(D33-X2*D54)
      FJ2= F1/R3-F3*D12+F3*X2*Y2*D54
      FJ3= A3/R3-FJ2
      FJ4=-F3*XY/R5-FJ1
      FK1=-Y*(D32-X2*D53)
      FK2=-X*(D32-Y2*D53)
      FK3=-F3*X*D/R5-FK2
# C---#
#C======================================
#C=====  STRIKE-SLIP CONTRIBUTION  =====
#C======================================
      if (POT1 != F0):
        DU[0]  =-X2*QR  -ALP3*FI1*SD
        DU[1]  =-XY*QR  -ALP3*FI2*SD
        DU[2]  =-C*X*QR -ALP3*FI4*SD
        DU[3]  =-X*QR*(F1+A5) -ALP3*FJ1*SD
        DU[4]  =-Y*QR*A5      -ALP3*FJ2*SD
        DU[5]  =-C*QR*A5      -ALP3*FK1*SD
        DU[6]  =-F3*X2/R5*UY      -ALP3*FJ2*SD
        DU[7]  =-F3*XY/R5*UY-X*QR -ALP3*FJ4*SD
        DU[8]  =-F3*C*X/R5*UY     -ALP3*FK2*SD
        DU[9]  =-F3*X2/R5*UZ  +ALP3*FK1*SD
        DU[10] =-F3*XY/R5*UZ  +ALP3*FK2*SD
        DU[11] = F3*X/R5*(-C*UZ +ALP3*Y*SD)
        for I in range(12):
           U[I]=U[I]+POT1/PI2*DU[I]
#      ENDIF
#C===================================
#C=====  DIP-SLIP CONTRIBUTION  =====
#C===================================
      if (POT2 != F0):
        DU[0]  = -X*P*QR +ALP3*FI3*SDCD
        DU[1]  =-Y*P*QR +ALP3*FI1*SDCD
        DU[2]  =-C*P*QR +ALP3*FI5*SDCD
        DU[3]  =-P*QR*A5 +ALP3*FJ3*SDCD
        DU[4]  = Y*P*QRX +ALP3*FJ1*SDCD
        DU[5]  = C*P*QRX +ALP3*FK3*SDCD
        DU[6]  = -F3*X/R5*VY      +ALP3*FJ1*SDCD
        DU[7]  = -F3*Y/R5*VY-P*QR +ALP3*FJ2*SDCD
        DU[8]  = -F3*C/R5*VY      +ALP3*FK1*SDCD
        DU[9]  = -F3*X/R5*VZ -ALP3*FK3*SDCD
        DU[10] =-F3*Y/R5*VZ -ALP3*FK1*SDCD
        DU[11] =-F3*C/R5*VZ +ALP3*A3/R3*SDCD
        for I in range(12):
           U[I]=U[I]+POT2/PI2*DU[I]
#      ENDIF
#C========================================
#C=====  TENSILE-FAULT CONTRIBUTION  =====
#C========================================
      if (POT3 != F0):
        DU[0]  = X*Q*QR -ALP3*FI3*SDSD
        DU[1]  = Y*Q*QR -ALP3*FI1*SDSD
        DU[2]  = C*Q*QR -ALP3*FI5*SDSD
        DU[3]  = Q*QR*A5 -ALP3*FJ3*SDSD
        DU[4]  =-Y*Q*QRX -ALP3*FJ1*SDSD
        DU[5]  =-C*Q*QRX -ALP3*FK3*SDSD
        DU[6]  = X*QR*WY     -ALP3*FJ1*SDSD
        DU[7]  = QR*(Y*WY+Q) -ALP3*FJ2*SDSD
        DU[8]  = C*QR*WY     -ALP3*FK1*SDSD
        DU[9]  = X*QR*WZ +ALP3*FK3*SDSD
        DU[10] = Y*QR*WZ +ALP3*FK1*SDSD
        DU[11] = C*QR*WZ -ALP3*A3/R3*SDSD
        for I in range(12):
           U[I]=U[I]+POT3/PI2*DU[I]
#      ENDIF
#C=========================================
#C=====  INFLATE SOURCE CONTRIBUTION  =====
#C=========================================
      if (POT4 != F0):
        DU[0]  = ALP3*X/R3
        DU[1]  = ALP3*Y/R3
        DU[2]  = ALP3*D/R3
        DU[3]  = ALP3*A3/R3
        DU[4]  =-ALP3*F3*XY/R5
        DU[5]  =-ALP3*F3*X*D/R5
        DU[6]  = DU[4]
        DU[7]  = ALP3*B3/R3
        DU[8]  =-ALP3*F3*Y*D/R5
        DU[9]  =-DU[5]
        DU[10] =-DU[8]
        DU[11] =-ALP3*C3/R3
        for I in range(12):
           U[I]=U[I]+POT4/PI2*DU[I]
      return U

def UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,\
        ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D, \
             P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3):
#      SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)
#      IMPLICIT REAL*8 (A-H,O-Z)
#      DIMENSION U(12),DU(12)
#C
#C********************************************************************
#C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
#C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
#C********************************************************************
#C
#C***** INPUT
#C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM
#C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
#C***** OUTPUT
#C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
#C
#     COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
#      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3
#     DATA F0,F1,F2,F3,F5,F7,F10,F15
#     *        /0.D0,1.D0,2.D0,3.D0,5.D0,7.D0,10.D0,15.D0/
#      DATA PI2/6.283185307179586D0/
      F0 = 0.0
      F1 = 1.0
      F2 = 2.0
      F3 = 3.0
      F5 = 5.0
      F7 = 7.0
      F10 = 10.0
      F15 = 15.0
      PI2 = np.pi*2
      #
      U = [np.copy(X)*F0]*12
      DU = [np.copy(X)*F0]*12
      C=D+Z
      #
      Q2=Q*Q
      R7=R5*R2
      A7=F1-F7*X2/R2
      B5=F1-F5*Y2/R2
      B7=F1-F7*Y2/R2
      C5=F1-F5*D2/R2
      C7=F1-F7*D2/R2
      D7=F2-F7*Q2/R2
      QR5=F5*Q/R2
      QR7=F7*Q/R2
      DR5=F5*D/R2
#C-----
      for I in range(12):
        U[I]=F0
# C======================================
# C=====  STRIKE-SLIP CONTRIBUTION  =====
# C======================================
      if (POT1 != F0):
        DU[0] =-ALP4*A3/R3*CD  +ALP5*C*QR*A5
        DU[1] = F3*X/R5*( ALP4*Y*CD +ALP5*C*(SD-Y*QR5) )
        DU[2] = F3*X/R5*(-ALP4*Y*SD +ALP5*C*(CD+D*QR5) )
        DU[3] = ALP4*F3*X/R5*(F2+A5)*CD   -ALP5*C*QRX*(F2+A7)
        DU[4] = F3/R5*( ALP4*Y*A5*CD +ALP5*C*(A5*SD-Y*QR5*A7) )
        DU[5] = F3/R5*(-ALP4*Y*A5*SD +ALP5*C*(A5*CD+D*QR5*A7) )
        DU[6] = DU[4]
        DU[7] = F3*X/R5*( ALP4*B5*CD -ALP5*F5*C/R2*(F2*Y*SD+Q*B7) )
        DU[8] = F3*X/R5*(-ALP4*B5*SD +ALP5*F5*C/R2*(D*B7*SD-Y*C7*CD) )
        DU[9] = F3/R5*   (-ALP4*D*A5*CD +ALP5*C*(A5*CD+D*QR5*A7) )
        DU[10]= F15*X/R7*( ALP4*Y*D*CD  +ALP5*C*(D*B7*SD-Y*C7*CD) )
        DU[11]= F15*X/R7*(-ALP4*Y*D*SD  +ALP5*C*(F2*D*CD-Q*C7) )
        #
        for I in range(12):
          U[I]=U[I]+POT1/PI2*DU[I]
#      ENDIF
#C===================================
#C=====  DIP-SLIP CONTRIBUTION  =====
#C===================================
      if (POT2 != F0):
        DU[0]  = ALP4*F3*X*T/R5          -ALP5*C*P*QRX
        DU[1]  =-ALP4/R3*(C2D-F3*Y*T/R2) +ALP5*F3*C/R5*(S-Y*P*QR5)
        DU[2]  =-ALP4*A3/R3*SDCD         +ALP5*F3*C/R5*(T+D*P*QR5)
        DU[3]  = ALP4*F3*T/R5*A5         -ALP5*F5*C*P*QR/R2*A7
        DU[4]  = F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-ALP5*F5*C/R2*(S-Y*P*QR7))
        DU[5]  = F3*X/R5*(ALP4*(F2+A5)*SDCD   -ALP5*F5*C/R2*(T+D*P*QR7))
        DU[6]  = DU[4]
        DU[7]  = F3/R5*(ALP4*(F2*Y*C2D+T*B5)\
     *                               +ALP5*C*(S2D-F10*Y*S/R2-P*QR5*B7))
        DU[8]  = F3/R5*(ALP4*Y*A5*SDCD-ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7))
        DU[9]  = F3*X/R5*(-ALP4*(S2D-T*DR5) -ALP5*F5*C/R2*(T+D*P*QR7))
        DU[10] = F3/R5*(-ALP4*(D*B5*C2D+Y*C5*S2D)\
     *                                -ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7))
        DU[11] = F3/R5*(-ALP4*D*A5*SDCD-ALP5*C*(S2D-F10*D*T/R2+P*QR5*C7))
        #
        for I in range(12):
          U[I]=U[I]+POT2/PI2*DU[I]
#       ENDIF
# C========================================
# C=====  TENSILE-FAULT CONTRIBUTION  =====
# C========================================
      if (POT3 != F0):
        DU[0] = F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))
        DU[1] = ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)
        DU[2] =-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)
        DU[3] =-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)
        DU[4] = F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)
     *                               -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z))
        DU[5] = F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)
     *                               +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z))
        DU[6] = DU[4]
        DU[7] = F3/R5*(-ALP4*(F2*Y*S2D+S*B5)
     *                -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5))
        DU[8] = F3/R5*( ALP4*Y*(F1-A5*SDSD)
     *                +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))
        DU[9] = F3*X/R5*(-ALP4*(C2D+S*DR5)
     *               +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))
        DU[10]= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)
     *               +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))
        DU[11]= F3/R5*(-ALP4*D*(F1-A5*SDSD)
     *               -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5)))
        for I in range(12):
           U[I]=U[I]+POT3/PI2*DU[I]
#       ENDIF
# C=========================================
# C=====  INFLATE SOURCE CONTRIBUTION  =====
# C=========================================
      if (POT4 != F0):
        DU[0]  = ALP4*F3*X*D/R5
        DU[1]  = ALP4*F3*Y*D/R5
        DU[2]  = ALP4*C3/R3
        DU[3]  = ALP4*F3*D/R5*A5
        DU[4]  = -ALP4*F15*XY*D/R7
        DU[5]  = -ALP4*F3*X/R5*C5
        DU[6]  = DU[4]
        DU[7]  = ALP4*F3*D/R5*B5
        DU[8]  = -ALP4*F3*Y/R5*C5
        DU[9]  = DU[5]
        DU[10] = DU[8]
        DU[11] = ALP4*F3*D/R5*(F2+C5)
        #
        for I in range(12):
          U[I]=U[I]+POT4/PI2*DU[I]
#       ENDIF
#       RETURN
#       ENDENDIF
# C========================================
# C=====  TENSILE-FAULT CONTRIBUTION  =====
# C========================================
      if (POT3 != F0):
        DU[0] = F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))
        DU[1] = ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)
        DU[2] =-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)
        DU[3] =-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)
        DU[4] = F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)\
     *                               -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z))
        DU[5] = F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)\
     *                               +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z))
        DU[6] = DU[4]
        DU[7] = F3/R5*(-ALP4*(F2*Y*S2D+S*B5)\
     *                -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5))
        DU[8] = F3/R5*( ALP4*Y*(F1-A5*SDSD)\
     *                +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))
        DU[9] = F3*X/R5*(-ALP4*(C2D+S*DR5)\
     *               +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))
        DU[10]= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)\
     *               +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))
        DU[11]= F3/R5*(-ALP4*D*(F1-A5*SDSD)\
     *               -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5)))
        #
        for I in range(12):
           U[I] = U[I]+POT3/PI2*DU[I]
#       ENDIF
# C=========================================
# C=====  INFLATE SOURCE CONTRIBUTION  =====
# C=========================================
      if (POT4 != F0):
        DU[0]  = ALP4*F3*X*D/R5
        DU[1]  = ALP4*F3*Y*D/R5
        DU[2]  = ALP4*C3/R3
        DU[3]  = ALP4*F3*D/R5*A5
        DU[4]  =-ALP4*F15*XY*D/R7
        DU[5]  =-ALP4*F3*X/R5*C5
        DU[6]  = DU[4]
        DU[7]  = ALP4*F3*D/R5*B5
        DU[8]  =-ALP4*F3*Y/R5*C5
        DU[9]  = DU[5]
        DU[10] = DU[8]
        DU[11] = ALP4*F3*D/R5*(F2+C5)
        #
        for I in range(12):
           U[I]=U[I]+POT4/PI2*DU[I]
#      ENDIF
      return U
#      END
###############################################################################
def d2rad(indegree):
    #
    return indegree/180*np.pi
###############################################################################    
def benchmark3d():
    x = np.array([0.5,0.5])
    y = np.array([2.6579,2.6579])
    z = np.array([0.,0.])
    U = okadaDC3D(0.6666667,x,y,z,3.060307,70.,-1.5,1.5,-1,1,1.,0.,0.)
    #
    print("  Output from strike-slip including EW: %-10.5f NS: %-10.5f and U: %-10.5f " % (U[0][0],\
                                     U[1][0],U[2][0]))
    U = okadaDC3D(0.6666667,x,y,z,3.060307,70.,-1.5,1.5,-1,1,0.,1.,0.)
    print("  Output from dip-slip    including EW: %-10.5f NS: %-10.5f and U: %-10.5f " % (U[0][0],\
                                     U[1][0],U[2][0]))     
    U = okadaDC3D(0.6666667,x,y,z,3.060307,70.,-1.5,1.5,-1,1,0.,0.,1.)
    print("  Output from tensile     including EW: %-10.5f NS: %-10.5f and U: %-10.5f " % (U[0][0],\
                                     U[1][0],U[2][0]))                                      
    return
############################################################################### 
# @numba.jit(parallel=True, nogil=True)   
#
def okadaDC3D(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3):
      '''
      DC3D, a fortran function was provided along with Okada's 1992 BSSA paper 
            a three-dimensional dislocation calculation. 
      This is the first python version that was transformed by Wanpeng Feng, @Ottawa, 
      on 2nd of October, 2016. The code was purely copied and modified based on the 
      original fortran codes from http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html
      #
      Note that the code has been improved to be used to support matrix computation.
      The NX1 coordinates vectors can be input simulatenously for surface deformation calculation.
      #
      All inputs are completely identical with the original fortran codes.
      #
      Warning: the function was only validated with my own matlab okada3D codes. 
               the return code (IRET) was temporarily out of use. This may be reconsidered again 
               at some of points in furture. Wanpeng Feng@Ottawa, Canada, 2016-10-03
               wanpeng.feng@hotmail.com
      
      '''
      '''
      SUBROUTINE  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,
      *              AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,
      *              UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4   ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,
      *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
      C
      C********************************************************************
      C*****                                                          *****
      C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****  
      C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
      C*****              CODED BY  Y.OKADA ... SEP.1991              *****
      C*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****
      C*****                          JUL.1993, MAY.2002              *****
      C********************************************************************
      C
      C***** INPUT
      C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
      C*****   X,Y,Z : COORDINATE OF OBSERVING POINT
      C*****   DEPTH : DEPTH OF REFERENCE POINT
      C*****   DIP   : DIP-ANGLE (DEGREE)
      C*****   AL1,AL2   : FAULT LENGTH RANGE
      C*****   AW1,AW2   : FAULT WIDTH RANGE
      C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
      C
      C***** OUTPUT
      C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)
      C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /
      C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )
      C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE
      C*****   IRET        : RETURN CODE
      C*****               :   =0....NORMAL
      C*****               :   =1....SINGULAR
      C*****               :   =2....POSITIVE Z WAS GIVEN
      C
      COMMON /C0/DUMMY(5),SD,CD
      DIMENSION  XI(2),ET(2),KXI(2),KET(2)
      DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)
      DATA  F0,EPS/ 0.D0, 1.D-6 /
      C----
      UX=U(1)
      UY=U(2)
      UZ=U(3)
      UXX=U(4)
      UYX=U(5)
      UZX=U(6)
      UXY=U(7)
      UYY=U(8)
      UZY=U(9)
      UXZ=U(10)
      UYZ=U(11)
      UZZ=U(12)-
      
      '''
      F0 = 0.
      EPS= 1.e-6
      #
      #IRET=0
      U   = [None] * 12
      if (np.sum(np.where(Z>0)) > 0.): 
         print(" WARNINNG: some of Z is over 0. Check again...")
         return
      #C-----
      #
      U      = [np.copy(X)*0.] * 12#[None] * 12
      DU     = [np.copy(X)*0.] * 12
      DUA    = [np.copy(X)*0.] * 12
      AALPHA = ALPHA
      DDIP   = DIP
      ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D = DCCON0(AALPHA,DDIP)
      # print(ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D)
      #C-----
      ZZ  = Z
      DD1 = DISL1
      DD2 = DISL2
      DD3 = DISL3
      #XI  = [X-AL1,X-AL2]
      XI0 = X-AL1
      # 
      XI0[np.abs(XI0)<EPS] = F0
      XI1 = X-AL2
      XI1[np.abs(XI1)<EPS] = F0 
      #
      XI = [XI0,XI1]
      #
      # C======================================
      # C=====  REAL-SOURCE CONTRIBUTION  =====
      # C======================================
      #
      D  = DEPTH+Z
      P  = Y*CD+D*SD
      Q  = Y*SD-D*CD
      #
      Q[np.abs(Q) < EPS] = F0
      #
      ET = [P-AW1,P-AW2]
      
      #if (np.abs(ET[0]) < EPS):
      ET[0][np.abs(ET[0])<EPS] = F0
      #if (np.abs(ET[1]) < EPS):
      ET[1][np.abs(ET[1])<EPS] = F0
      #C--------------------------------
      #C----- REJECT SINGULAR CASE -----
      #C--------------------------------
      #C----- ON FAULT EDGE
      #if (Q == F0 and ((XI[0]*XI[1] < F0 and  ET[0]*ET[1] == F0) or \
      #   (ET[0]*ET[1] <= F0 and  XI[0]*XI[1] == F0))):
      # IRET=1
      #  return U*0.
      #
      # print('TEST!!! %f %f %f %f ' % (XI[0],XI[1],ET[0],ET[1]))
      # C----- ON NEGATIVE EXTENSION OF FAULT EDGE
      KXI = [XI[0]*0.,XI[0]*0]
      #KXI(2)=0
      KET = [XI[0]*0.,XI[0]*0]
      #KET(2)=0
      R12 = np.sqrt(XI[0]*XI[0]+ET[1]*ET[1]+Q*Q)
      R21 = np.sqrt(XI[1]*XI[1]+ET[0]*ET[0]+Q*Q)
      R22 = np.sqrt(XI[1]*XI[1]+ET[1]*ET[1]+Q*Q)
      #
      #if (XI[0] < F0 and R21+XI[1] < EPS):
      # print(R21.shape)
      KXI[0][(XI[0] < F0) & (R21+XI[1] < EPS)] = 1.
      #if (XI[0]< F0 and R22+XI[1] < EPS):
      KXI[1][(XI[0] < F0) & (R22+XI[1] < EPS)] = 1.
      #if (ET[0] < F0 and R12+ET[1] < EPS):
      KET[0][(ET[0] < F0) & (R12+ET[1] < EPS)] = 1.
      #if (ET[0] < F0 and R22+ET[1] < EPS):
      KET[1][(ET[0] < F0) & (R22+ET[1] < EPS)] = 1.
#C=====
      #DO 223 K=1,2
      #DO 222 J=1,2
      counter = 0
      for K in range(2):
          for J in range(2):
             # 
             counter += 1
             XI2,ET2,Q2,R,R2,R3,R5,TEMPY,TEMPD,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ\
                 = DCCON2(XI[J],ET[K],Q,SD,CD,KXI[K],KET[J])
             #           
             DUA = UA(XI[J],ET[K],Q,DD1,DD2,DD3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                     XI2,ET2,Q2,R,R2,R3,R5,TEMPY,TEMPD,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ)
             #
             #
             for I in range(0,10,3):
               #
               DU[I]   =  DUA[I] * -1.
               DU[I+1] = -1.*DUA[I+1]*CD+DUA[I+2]*SD
               DU[I+2] = -1.*DUA[I+1]*SD-DUA[I+2]*CD
               if (I >= 9.):
                 # print(" TeST: %f %f %f %f" % (DU[10],DUA[10],CD,SD)) 
                 DU[I]   = -1.*DU[I]
                 DU[I+1] = -1.*DU[I+1] 
                 DU[I+2] = -1.*DU[I+2]
               
     #  220   CONTINUE
               
             for I in range(12):
                if (J+K != 1):
                   if counter == 1:
                      U[I] = DU[I]*1.0
                   else:
                      U[I] = U[I]+DU[I]
                if (J+K == 1):
                    if counter == 1:
                      U[I] = DU[I]*1.0
                    else:
                      U[I] = U[I]-DU[I]
      #C=======================================
      #C=====  IMAGE-SOURCE CONTRIBUTION  =====
      #C=======================================
      
      D = DEPTH-Z
      P = Y*CD+D*SD
      Q = Y*SD-D*CD
      # 
      ET = [P-AW1,P-AW2]
      #
      #if (np.abs(Q) < EPS):
      Q[np.abs(Q)<EPS] = F0
      #if (np.abs(ET[0]) < EPS):
      ET[0][np.abs(ET[0]) < EPS] = F0
      #if (np.abs(ET[1]) < EPS):
      ET[1][np.abs(ET[0]) < EPS] = F0
      #C--------------------------------
      #C----- REJECT SINGULAR CASE -----
      #C--------------------------------
      #C----- ON FAULT EDGE
      KXI=[X*0.,X*0.]
      #KXI(2)=0
      KET=[X*0.,X*0.]  
      #
      R12=np.sqrt(XI[0]*XI[0]+ET[1]*ET[1]+Q*Q)
      R21=np.sqrt(XI[1]*XI[1]+ET[0]*ET[0]+Q*Q)
      R22=np.sqrt(XI[1]*XI[1]+ET[1]*ET[1]+Q*Q)
      #if (XI[0] < F0 and R21+XI[1] < EPS):
      KXI[0][(XI[0] < F0) & (R21+XI[1] < EPS)] = 1.
      #if (XI[0] < F0 and R22+XI[1] < EPS):
      KXI[1][(XI[0] < F0) & (R22+XI[1] < EPS)] = 1.
      #if (ET[0] < F0 and R12+ET[1] < EPS):
      KET[0][(ET[0] < F0) & (R12+ET[1] < EPS)] = 1.
      #if (ET[0] < F0 and R22+ET[1] < EPS):
      KET[1][(ET[0] < F0) & (R22+ET[1] < EPS)] = 1
      #C=====
      for K in range(2):
         for J in range(2):
             XI2,ET2,Q2,R,R2,R3,R5,TY,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ\
                 = DCCON2(XI[J],ET[K],Q,SD,CD,KXI[K],KET[J])
             #print(Q,R)      
             DUA = UA(XI[J],ET[K],Q,DD1,DD2,DD3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                      XI2,ET2,Q2,R,R2,R3,R5,TY,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ)
             DUB = UB(XI[J],ET[K],Q,DD1,DD2,DD3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                      XI2,ET2,Q2,R,R2,R3,R5,TY,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ)
             #
             DUC = UC(XI[J],ET[K],Q,ZZ,DD1,DD2,DD3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
                      XI2,ET2,Q2,R,R2,R3,R5,TY,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ)
             # print(DUC[-6:])       
             #C-----
             for I in range(0,10,3):
                # 
                DU[I]   = DUA[I]+DUB[I]+Z*DUC[I]
                DU[I+1] = (DUA[I+1]+DUB[I+1]+Z*DUC[I+1])*CD-(DUA[I+2]+DUB[I+2]+Z*DUC[I+2])*SD
                DU[I+2] = (DUA[I+1]+DUB[I+1]-Z*DUC[I+1])*SD+(DUA[I+2]+DUB[I+2]-Z*DUC[I+2])*CD
                if (I >= 9):
                  DU[ I]  = DU[I] + DUC[0]
                  DU[I+1] = DU[I+1]+DUC[1]*CD-DUC[2]*SD
                  DU[I+2] = DU[I+2]-DUC[1]*SD-DUC[2]*CD
             #  330   CONTINUE
             # print(DU[0:6])
             #for I in range(12):
             for I in range(12):
                if (J+K != 1):
                   U[I]=U[I]+DU[I]
                if (J+K == 1):
                   U[I]=U[I]-DU[I]
    
      return U
      # the last line of DC3D
###############################################################################
# @numba.jit(parallel=True, nogil=True)      
def UA(XI,ET,Q,DISL1,DISL2,DISL3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
       XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ):
      '''
      SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(12),DU(12)
      C
      C********************************************************************
      C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****
      C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
      C********************************************************************
      C
      C***** INPUT
      C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
      C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
      C***** OUTPUT
      C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
      C
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,
      *           EY,EZ,FY,FZ,GY,GZ,HY,HZ
      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/
      C-----
      '''
      F0 = 0.
      F2 = 2.0
      PI2= np.pi * 2
      #
      U = [np.copy(XI) * 0] * 12
      DU = [np.copy(XI) * 0] * 12
      #
      XY = XI*Y11
      QX = Q *X11
      QY = Q *Y11
#C======================================
#C=====  STRIKE-SLIP CONTRIBUTION  =====
#C======================================
      if (DISL1 != F0):
        DU[ 0] =    TT/F2 +ALP2*XI*QY
        DU[ 1] =           ALP2*Q/R
        DU[ 2] = ALP1*ALE -ALP2*Q*QY
        DU[ 3] =-ALP1*QY  -ALP2*XI2*Q*Y32
        DU[ 4] =          -ALP2*XI*Q/R3
        DU[ 5] = ALP1*XY    + ALP2*XI*Q2*Y32
        DU[ 6] = ALP1*XY*SD + ALP2*XI*FY+D/F2*X11
        DU[ 7] =                    ALP2*EY
        DU[ 8] = ALP1*(CD/R+QY*SD) -ALP2*Q*FY
        DU[ 9] = ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11
        DU[10] =                    ALP2*EZ
        DU[11] =-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ
        for I in range(12):
            U[I]=DISL1/PI2*DU[I]
            #
#C======================================
#C=====    DIP-SLIP CONTRIBUTION   =====
#C======================================
      if (DISL2 != F0):
        DU[ 0] =           ALP2*Q/R
        DU[ 1] =    TT/F2 +ALP2*ET*QX
        DU[ 2] = ALP1*ALX -ALP2*Q*QX
        DU[ 3] =        -ALP2*XI*Q/R3
        DU[ 4] = -QY/F2 -ALP2*ET*Q/R3
        DU[ 5] = ALP1/R +ALP2*Q2/R3
        DU[ 6] =                      ALP2*EY
        DU[ 7] = ALP1*D*X11+XY/F2*SD +ALP2*ET*GY
        DU[ 8] = ALP1*Y*X11          -ALP2*Q*GY
        DU[ 9] =                      ALP2*EZ
        DU[10] = ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ
        DU[11]=-ALP1*D*X11          -ALP2*Q*GZ
        for I in range(12):
            if DISL1 != F0:
               U[I]=U[I]+DISL2/PI2*DU[I]
            else:
               U[I]=DISL2/PI2*DU[I]
        #    
#C========================================
#C=====  TENSILE-FAULT CONTRIBUTION  =====
#C========================================
      if (DISL3 != F0): 
        DU[ 0] =-ALP1*ALE -ALP2*Q*QY
        DU[ 1] =-ALP1*ALX -ALP2*Q*QX
        DU[ 2] =    TT/F2 -ALP2*(ET*QX+XI*QY)
        DU[ 3] =-ALP1*XY  +ALP2*XI*Q2*Y32
        DU[ 4] =-ALP1/R   +ALP2*Q2/R3
        DU[ 5] =-ALP1*QY  -ALP2*Q*Q2*Y32
        DU[ 6] =-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY
        DU[ 7] =-ALP1*Y*X11         -ALP2*Q*GY
        DU[ 8] = ALP1*(D*X11+XY*SD) +ALP2*Q*HY
        DU[ 9] = ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ
        DU[10] = ALP1*D*X11         -ALP2*Q*GZ
        DU[11] = ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ
        for I in range(12):
            if (DISL1 != F0 or DISL2 != F0):
               U[I]=U[I]+DISL3/PI2*DU[I]
            else:
               U[I]=DISL3/PI2*DU[I]
      return U
###############################################################################
# @numba.jit(parallel=True, nogil=True)      
def UB(XI,ET,Q,DISL1,DISL2,DISL3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
       XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ):
      '''
      SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(12),DU(12)
C
C********************************************************************
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
C********************************************************************
C
C***** INPUT
C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
C***** OUTPUT
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
C
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ
      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/
C-----
      '''
      F0  = 0.
      F1  = 1.
      F2  = 2.0
      PI2 = np.pi * 2
      RD  = R+D
      D11 = F1/(R*RD)
      AJ2 = XI*Y/RD*D11
      AJ5 = -(D+Y*Y/RD)*D11
      if (CD != F0):
        AI4 = XI*0
        AI4[XI==F0] = F0
        X           = np.sqrt(XI2+Q2)
        AI4[XI!=F0] = F1/CDCD*( XI[XI!=F0] / RD[XI!=F0]*SDCD +\
                        F2 * np.arctan((ET[XI!=F0]*(X[XI!=F0]+Q[XI!=F0]*CD)+\
                        X[XI!=F0]*(R[XI!=F0]+X[XI!=F0])*SD)/(XI[XI!=F0]*(R[XI!=F0]+X[XI!=F0])*CD)))
        #if (XI == F0):
        #  AI4=F0
        #else:
        #  X   = np.sqrt(XI2+Q2)
        #  AI4 = F1/CDCD*( XI/RD*SDCD+F2*np.arctan((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )
        #
        AI3   = (Y*CD/RD-ALE+SD*np.log(RD))/CDCD
        AK1   = XI*(D11-Y11*SD)/CD
        AK3   = (Q*Y11-Y*D11)/CD
        AJ3   = (AK1-AJ2*SD)/CD
        AJ6   = (AK3-AJ5*SD)/CD
      else:
        RD2=RD*RD
        AI3=(ET/RD+Y*Q/RD2-ALE)/F2
        AI4=XI*Y/RD2/F2
        AK1=XI*Q/RD*D11
        AK3=SD/RD*(XI2*D11-F1)
        AJ3=-XI/RD2*(Q2*D11-F1/F2)
        AJ6=-Y/RD2*(XI2*D11-F1/F2)
      #
#C-----
      XY  =  XI*Y11
      AI1 = -XI/RD*CD-AI4*SD
      AI2 =  np.log(RD)+AI3*SD
      AK2 =  F1/R+AK3*SD
      AK4 =  XY*CD-AK1*SD
      AJ1 =  AJ5*CD-AJ6*SD
      AJ4 = -XY-AJ2*CD+AJ3*SD
#C=====
      U = [np.copy(X11)*0] * 12
      DU = [np.copy(X11)*0] * 12
      #
      QX = Q*X11
      QY = Q*Y11
#C======================================
#C=====  STRIKE-SLIP CONTRIBUTION  =====
#C======================================
      if (DISL1 != F0): 
        DU[ 0] = -XI*QY-TT -ALP3*AI1*SD
        DU[ 1] = -Q/R      +ALP3*Y/RD*SD
        DU[ 2] =  Q*QY     -ALP3*AI2*SD
        DU[ 3] = XI2*Q*Y32 -ALP3*AJ1*SD
        DU[ 4] = XI*Q/R3   -ALP3*AJ2*SD
        DU[ 5] =-XI*Q2*Y32 -ALP3*AJ3*SD
        DU[ 6] =-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD
        DU[ 7] =-EY          +ALP3*(F1/R+AJ5)*SD
        DU[ 8] = Q*FY        -ALP3*(QY-AJ6)*SD
        DU[ 9] =-XI*FZ-Y*X11 +ALP3*AK1*SD
        DU[10] =-EZ          +ALP3*Y*D11*SD
        DU[11] = Q*FZ        +ALP3*AK2*SD
        #
        for I in range(12):
            U[I] = DISL1/PI2*DU[I]
            
#C======================================
#C=====    DIP-SLIP CONTRIBUTION   =====
#C======================================
      if (DISL2 != F0): 
        DU[ 0] =-Q/R         + ALP3*AI3*SDCD
        DU[ 1] =-ET*QX-TT    - ALP3*XI/RD*SDCD
        DU[ 2] = Q*QX        + ALP3*AI4*SDCD
        DU[ 3] = XI*Q/R3     + ALP3*AJ4*SDCD
        DU[ 4] = ET*Q/R3+QY  + ALP3*AJ5*SDCD
        DU[ 5] =-Q2/R3       + ALP3*AJ6*SDCD
        DU[ 6] =-EY          + ALP3*AJ1*SDCD
        DU[ 7] =-ET*GY-XY*SD + ALP3*AJ2*SDCD
        DU[ 8] = Q*GY        + ALP3*AJ3*SDCD
        DU[ 9] =-EZ          - ALP3*AK3*SDCD
        DU[10] =-ET*GZ-XY*CD - ALP3*XI*D11*SDCD
        DU[11] = Q*GZ        - ALP3*AK4*SDCD
        for I in range(12):
            if DISL1 != F0:
               U[I]=U[I]+DISL2/PI2*DU[I]
            else:
               U[I]=DISL2/PI2*DU[I]
            
#C========================================
#C=====  TENSILE-FAULT CONTRIBUTION  =====
#C========================================
      if (DISL3 != F0):
        DU[ 0] = Q*QY           -ALP3*AI3*SDSD
        DU[ 1] = Q*QX           +ALP3*XI/RD*SDSD
        DU[ 2] = ET*QX+XI*QY-TT -ALP3*AI4*SDSD
        DU[ 3] =-XI*Q2*Y32 -ALP3*AJ4*SDSD
        DU[ 4] =-Q2/R3     -ALP3*AJ5*SDSD
        DU[ 5] = Q*Q2*Y32  -ALP3*AJ6*SDSD
        DU[ 6] = Q*FY -ALP3*AJ1*SDSD
        DU[ 7] = Q*GY -ALP3*AJ2*SDSD
        DU[ 8] =-Q*HY -ALP3*AJ3*SDSD
        DU[ 9] = Q*FZ +ALP3*AK3*SDSD
        DU[10] = Q*GZ +ALP3*XI*D11*SDSD
        DU[11] =-Q*HZ +ALP3*AK4*SDSD
        for I in range(12):
            if (DISL1 != F0 or DISL2 != F0):
               U[I]=U[I]+DISL3/PI2*DU[I]
            else:
               U[I]=DISL3/PI2*DU[I]
            
      return U
#      END

###############################################################################
# @numba.jit(nopython=True, parallel=True, nogil=True)     
def UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D,\
       XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ):
      '''
      SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U,)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(12),DU(12)
C
C********************************************************************
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****
C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
C********************************************************************
C
C***** INPUT
C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM
C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
C***** OUTPUT
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
C
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,
     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ
      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/
C-----
      '''
      F0 = 0.
      F1 = 1.
      F2 = 2.0
      F3 = 3.0
      PI2 = np.pi * 2
      C   = D+Z
      X53 = (8.0*R2+9.0*R*XI+F3*XI2)*X11*X11*X11/R2
      Y53 = (8.0*R2+9.0*R*ET+F3*ET2)*Y11*Y11*Y11/R2
      H   = Q*CD-Z
      Z32 = SD/R3-H*Y32
      Z53 = F3*SD/R5-H*Y53
      Y0  = Y11-XI2*Y32
      Z0  = Z32-XI2*Z53
      PPY = CD/R3+Q*Y32*SD
      PPZ = SD/R3-Q*Y32*CD
      QQ  = Z*Y32+Z32+Z0
      QQY = F3*C*D/R5-QQ*SD
      QQZ = F3*C*Y/R5-QQ*CD+Q*Y32
      XY  = XI*Y11
      # QX=Q*X11
      QY  = Q*Y11
      QR  = F3*Q/R5
      # CQX=C*Q*X53
      CDR = (C+D)/R3
      YY0 = Y/R3-Y0*CD
#C=====
#      DO 111  I=1,12
#  111 U(I)=F0
      U    = [np.copy(XI)*0.] * 12
      DU   = [np.copy(XI)*0.] * 12
#C======================================
#C=====  STRIKE-SLIP CONTRIBUTION  =====
#C======================================
      if (DISL1 != F0): 
        DU[ 0] = ALP4*XY*CD           -ALP5*XI*Q*Z32
        DU[ 1] = ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3
        DU[ 2] = ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)
        DU[ 3] = ALP4*Y0*CD                  -ALP5*Q*Z0
        DU[ 4] =-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR
        DU[ 5] =-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)
        DU[ 6] =-ALP4*XI*PPY*CD    -ALP5*XI*QQY
        DU[ 7] = ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD-ALP5*(CDR*SD-ET/R3-C*Y*QR)
        DU[ 8] =-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD)
        DU[ 9] = ALP4*XI*PPZ*CD    -ALP5*XI*QQZ
        DU[10] = ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR)
        DU[11] = YY0*C - ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)
        #for I in range(12):
        for I in range(12):
            U[I]=DISL1/PI2*DU[I]
#C======================================
#C=====    DIP-SLIP CONTRIBUTION   =====
#C======================================
      if (DISL2 != F0):
        DU[ 0] = ALP4*CD/R -QY*SD -ALP5*C*Q/R3
        DU[ 1] = ALP4*Y*X11       -ALP5*C*ET*Q*X32
        DU[ 2] =     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32)
        DU[ 3] =-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD
        DU[ 4] =-ALP4*Y/R3     +ALP5*C*ET*QR
        DU[ 5] =    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2)
        DU[ 6] =-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR)
        DU[ 7] = ALP4*(X11-Y*Y*X32) -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)
        DU[ 8] =  XI*PPY*SD+Y*D*X32 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)
        DU[ 9] =      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR)
        DU[10] = ALP4*Y*D*X32       -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)
        DU[11] =-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)
        for I in range(12):
            if DISL1 != F0:
               U[I]=U[I]+DISL2/PI2*DU[I]
            else:
               U[I]=DISL2/PI2*DU[I]
        
#C========================================
#C=====  TENSILE-FAULT CONTRIBUTION  =====
#C========================================
      if (DISL3 != F0):
        DU[ 0] = -ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32)
        DU[ 1] = ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32)
        DU[ 2] = ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32)
        DU[ 3] = ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)
        DU[ 4] = ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2)
        DU[ 5] =-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0)
        DU[ 6] = ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)
        DU[ 7] =-ALP4*F2*XI*PPY*SD-Y*D*X32+ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)
        DU[ 8] =-ALP4*(XI*PPY*CD-X11+Y*Y*X32)+ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY)
        DU[ 9] =  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)
        DU[10] = ALP4*F2*XI*PPZ*SD-X11+D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)
        DU[11] = ALP4*(XI*PPZ*CD+Y*D*X32)+ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ)
        for I in range(12):
            if (DISL1 != F0 or DISL2 != F0):
               U[I]=U[I]+DISL3/PI2*DU[I]
            else:
               U[I]=DISL3/PI2*DU[I]
            
      return U
#      END
###############################################################################
# @numba.jit(nopython=True, parallel=True, nogil=True)      
def DCCON0(ALPHA,DIP):
      '''
      SUBROUTINE  DCCON0(ALPHA,DIP)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C*******************************************************************
C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****
C*******************************************************************
C
C***** INPUT
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
C*****   DIP   : DIP-ANGLE (DEGREE)
C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO
C
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
      DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/
      DATA EPS/1.D-6/
      '''
      F0 = 0.
      F1 = 1.e+0
      F2 = 2.e+0
      PI2= np.pi*2
      EPS= 1.e-6
      
#C-----
      ALP1 = (F1-ALPHA)/F2
      ALP2 = ALPHA/F2
      ALP3 = (F1-ALPHA)/ALPHA
      ALP4 = F1-ALPHA
      ALP5 = ALPHA
#C-----
      P18  = PI2/360.0
      SD   = np.sin(DIP*P18)
      CD   = np.cos(DIP*P18)
      if (np.abs(CD) < EPS):
        CD = F0
        if (SD > F0):
            SD= F1
        if (SD < F0):
            SD=-F1
      #
      SDSD = SD*SD
      CDCD = CD*CD
      SDCD = SD*CD
      S2D  = F2*SDCD
      C2D  = CDCD-SDSD
      return ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D
###############################################################################
# @numba.jit(nopython=True, parallel=True, nogil=True)      
def DCCON1(X,Y,D,SD,CD):
      '''      
      SUBROUTINE  DCCON1(X,Y,D)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**********************************************************************
C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    *****
C**********************************************************************
C
C***** INPUT
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM
C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO
C
      COMMON /C0/DUMMY(5),SD,CD
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,
     *           UY,VY,WY,UZ,VZ,WZ
      DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/
      '''
      F0 = 0.
      F1 = 1.e+0
      F3 = 3.e+0
      F5 = 5.e+0
      EPS= 1.e-6
      
#C-----
      # if (np.abs(X) < EPS):
      X[np.abs(X)<EPS] = F0
      #if (np.abs(Y) < EPS):
      Y[np.abs(Y)<EPS] = F0
      #if (np.abs(D) < EPS):
      D[np.abs(D)<EPS] = F0
      P  = Y*CD+D*SD
      Q  = Y*SD-D*CD
      S  = P*SD+Q*CD
      T  = P*CD-Q*SD
      XY = X*Y
      X2 = X*X
      Y2 = Y*Y
      D2 = D*D
      R2 = X2+Y2+D2
      R  = np.sqrt(R2)
      if np.sum(np.where(R==F0)) != 0:
         print(" ERROR from DCCON1: some or all of R  is 0!!!")
         return
      R3 = R *R2
      R5 = R3*R2
      # R7=R5*R2
#C-----
      A3 = F1-F3*X2/R2
      A5 = F1-F5*X2/R2
      B3 = F1-F3*Y2/R2
      C3 = F1-F3*D2/R2
#C-----
      QR  = F3*Q/R5
      QRX = F5*QR*X/R2
#C-----
      UY  = SD-F5*Y*Q/R2
      UZ  = CD+F5*D*Q/R2
      VY  = S-F5*Y*P*Q/R2
      VZ  = T+F5*D*P*Q/R2
      WY  = UY+SD
      WZ  = UZ+CD
      return P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ
###############################################################################
# @numba.jit(nopython=True, parallel=True, nogil=True)     
#      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD,KXI,KET)
def DCCON2(XI,ET,Q,SD,CD,KXI,KET):
      '''
#      IMPLICIT REAL*8 (A-H,O-Z)
#C
#C**********************************************************************
#C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   *****
#C**********************************************************************
#C
#C***** INPUT
#C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
#C*****   SD,CD   : SIN, COS OF DIP-ANGLE
#C*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
#C
#C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0
#C
#      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,
#     *           EY,EZ,FY,FZ,GY,GZ,HY,HZ
#      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/
#C-----
      '''
      F0 = 0.
      F1 = 1.
      F2 = 2.
      EPS = 1.e-6
      #
      #if (np.abs(XI) < EPS):
      XI[np.abs(XI) < EPS] = F0
      #if (np.abs(ET) < EPS):
      ET[np.abs(ET) < EPS] = F0
      #if (np.abs(Q) < EPS):
      Q[np.abs(Q) < EPS] = F0
      #    
      XI2 = XI*XI
      ET2 = ET*ET
      Q2  = Q*Q
      R2  = XI2+ET2+Q2
      R   = np.sqrt(R2)
      #if (R == F0):
      #    return
      R3 = R *R2
      R5 = R3*R2
      Y  = ET*CD+Q*SD
      D  = ET*SD-Q*CD
#C-----
      #if (Q == F0):
      TT        = XI*0
      # TT=F0
      TT[Q==F0] = F0
      TT[Q!=F0] = np.arctan(XI[Q!=F0]*ET[Q!=F0]/(Q[Q!=0]*R[Q!=0]))
      # TT=np.arctan(XI*ET/(Q*R))
#C-----
      #if (KXI == 1):
      ALX         = KXI*0
      RXI         = R + XI
      ALX[KXI==1] = -np.log(R[KXI==1]-XI[KXI==1])
      ALX[KXI!=1] =  np.log(RXI[KXI!=1])
      X11         = KXI*0
      X11[KXI==1] = F0
      X11[KXI!=1] = F1 / (R[KXI!=1]*RXI[KXI!=1])
      X32         = KXI*0
      X32[KXI==1] = F0
      X32[KXI!=1] = (R[KXI!=1]+RXI[KXI!=1])*X11[KXI!=1]*X11[KXI!=1]/R[KXI!=1]
      #else:
      #  RXI=R+XI
      #  ALX=np.log(RXI)
      #  X11=F1/(R*RXI)
      #  X32=(R+RXI)*X11*X11/R
      #C-----
      #if (KET == 1):
      ALE         = KET*0
      RET         = R+ET
      ALE[KET==1] = -np.log(R[KET==1]-ET[KET==1])
      ALE[KET!=1] =  np.log(RET[KET!=1])
      Y11         = KET*0
      Y11[KET==1] = F0
      Y11[KET!=1] = F1/(R[KET!=1]*RET[KET!=1])
      Y32         = KET*0
      Y32[KET==1] = F0
      Y32[KET!=1] = (R[KET!=1]+RET[KET!=1])*(Y11[KET!=1]**2)/R[KET!=1]

#C-----
      EY=SD/R-Y*Q/R3
      EZ=CD/R+D*Q/R3
      FY=D/R3+XI2*Y32*SD
      FZ=Y/R3+XI2*Y32*CD
      GY=F2*X11*SD-Y*Q*X32
      GZ=F2*X11*CD+D*Q*X32
      HY=D*Q*X32+XI*Q*Y32*SD
      HZ=Y*Q*X32+XI*Q*Y32*CD
      return XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ
      # return XI2,ET2,Q2,R,R2,R3,R5,TT,ALX,ALE,Y11,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ,X11,X32

###############################################################################
def benchmark2d():
  #
  #
  alp         = 0.5
  x,y,dep,dip = 2.,3.,4.,70.
  #
  strtest  = okadaDC2D(alp,x,y,dep,0.,\
                     3.,0.,2.,np.sin(dip/180.*np.pi),\
                     np.cos(dip/180.*np.pi),1.,0.,0.)
                     
  print(" [1] Targeted results for Strike-Slip Model in Test1: -8.689e3  -4.298e-3  -2.747e-3 ")
  #
  print("     Calculated results for Test1:                    %-10f %-10f %-10f" % (strtest[0],strtest[1],strtest[2]))
  #
  strtest  = okadaDC2D(alp,x,y,dep,0.,\
                     3.,0.,2.,np.sin(dip/180.*np.pi),\
                     np.cos(dip/180.*np.pi),0.,1.,0.)
                     
  print(" [2] Targeted results for Dip-Slip Model in Test1: -4.682e-3  -3.527e-3  -3.564e-2 ")
  #
  print("     Calculated results for Test1:                    %-10f %-10f %-10f" % (strtest[0],strtest[1],strtest[2]))
  strtest  = okadaDC2D(alp,x,y,dep,0.,\
                     3.,0.,2.,np.sin(dip/180.*np.pi),\
                     np.cos(dip/180.*np.pi),0.,0.,1.)
                     
  print(" [3] Targeted results for Tensile-Slip Model in Test1: -2.660e-4  1.056e-2  3.214e-3 ")
  #
  print("     Calculated results for Test1:   %-10f %-10f %-10f" % (strtest[0],strtest[1],strtest[2]))
  return True
  #
###############################################################################  
def stresstensorrot(strike,stresstensor,ver=np.pi/2):
    #
    # To rotate a stress tensor from fault projection to the north-south and east-west
    # direction.
    #
    dstr = strike/180.*np.pi
    cosb = np.cos(dstr)
    sinb = np.sin(dstr)
    nn   = stresstensor.shape[0]
    #%
    t    = np.zeros([6,6])
    sn   = np.zeros([nn,6])
    #
    #
    bt = np.arcsin(sinb)
    #
    if cosb > 0.0:
       #% rotation of axes (positive:clock wise, negative, anti-cw)
       xbeta = -bt
       xdel  = 0.0
       ybeta = -bt + ver
       ydel  = 0.0
       zbeta = -bt - ver
       zdel  = ver
    else:
       xbeta = bt - np.pi
       xdel  = 0.0
       ybeta = bt - ver
       ydel  = 0.0;
       zbeta = bt - ver
       zdel  = ver  
    # end
    #
    xl = np.cos(xdel) * np.cos(xbeta)
    xm = np.cos(xdel) * np.sin(xbeta)
    xn = np.sin(xdel)

    yl = np.cos(ydel) * np.cos(ybeta)
    ym = np.cos(ydel) * np.sin(ybeta)
    yn = np.sin(ydel)
    #%
    zl = np.cos(zdel) * np.cos(zbeta)
    zm = np.cos(zdel) * np.sin(zbeta)
    zn = np.sin(zdel)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t[0,0] = xl * xl;
    t[0,1] = xm * xm;
    t[0,2] = xn * xn;
    t[0,3] = 2.0 * xm * xn
    t[0,4] = 2.0 * xn * xl
    t[0,5] = 2.0 * xl * xm
    t[1,0] = yl * yl
    t[1,1] = ym * ym
    t[1,3] = yn * yn
    t[1,3] = 2.0 * ym * yn
    t[1,4] = 2.0 * yn * yl
    t[1,5] = 2.0 * yl * ym
    t[2,0] = zl * zl
    t[2,1] = zm * zm
    t[2,2] = zn * zn
    t[2,3] = 2.0 * zm * zn
    t[2,4] = 2.0 * zn * zl
    t[2,5] = 2.0 * zl * zm
    t[3,0] = yl * zl
    t[3,1] = ym * zm
    t[3,2] = yn * zn  
    t[3,3] = ym * zn + zm * yn
    t[3,4] = yn * zl + zn * yl
    t[3,5] = yl * zm + zl * ym
    t[4,0] = zl * xl
    t[4,1] = zm * xm
    t[4,2] = zn * xn
    t[4,3] = xm * zn + zm * xn
    t[4,4] = xn * zl + zn * xl
    t[4,5] = xl * zm + zl * xm
    t[5,0] = xl * yl
    t[5,1] = xm * ym
    t[5,2] = xn * yn
    t[5,3] = xm * yn + ym * xn
    t[5,4] = xn * yl + yn * xl
    t[5,5] = xl * ym + yl * xm
    #
    for k in range(nn):
       sn[k,:] = np.dot(t[:,:],stresstensor[k,:])
    #  
    return sn
###############################################################################  
#
def fpara2stress(mfpara,x,y,z,alpha=0.6666666666666667,strainonly=False):
    #
    for ni in range(mfpara.shape[0]):
        if ni == 0:
            xx,yy,zz,xy,xz,yz = dc3dtostress(mfpara[ni,:],x,y,z,\
                                             alpha=alpha,\
                                             strainonly=strainonly)
        else:
            cxx,cyy,czz,cxy,cxz,cyz = dc3dtostress(mfpara[ni,:],x,y,z,\
                                                   alpha=alpha,\
                                                   strainonly=strainonly)
            xx = xx + cxx
            yy = yy + cyy
            zz = zz + czz
            xy = xy + cxy
            yz = yz + cyz
            xz = xz + cxz
        #
    return xx,yy,zz,xy,yz,xz
#
###############################################################################
# @numba.jit(nopython=True, parallel=True)    
def dc3dtostress(fpara,x,y,z,alpha=0.66666666666666667,young=1e6,pois=0.25,\
                 strainonly=False):
  '''
  strainonly can be applied to return strain tensor only...
  '''
  # Note that alpha used for okada3D is different with that for okada2D
  # 
  if not isinstance(x,np.ndarray):
  #if type(x) is not np.ndarray:
     #
     # python version dc3d supports matrix operation, which can largely make 
     # computation faster. The input coordinates, (x,y,z) all should be 
     # in numpy.ndarray. If it was found not numpy.array, 
     # here we will force them to be ...
     # by Wanpeng Feng, @Ottawa, 2016-10-02
     #
     x  = np.array([x])
     y  = np.array([y])
     z  = np.array([z])
  #
  fpara      = np.reshape(np.array(fpara),10)
  fpara[3]   = (fpara[3] == 90.) * 89.9999 + (fpara[3] != 90. ) * fpara[3]
  fx,fy      = geo2okada3d(x,y,fpara)
  #
  # % alpha is media parameter. 
  diprad     = fpara[3]/180.*np.pi
  DEP        = fpara[4] + fpara[5]/2 * np.sin(diprad)
  DIP        = fpara[3]
  strike     = fpara[2]
  AW1        = fpara[5]/2 * -1
  AW2        = fpara[5]/2
  AL1        = fpara[6]/2 * -1
  AL2        = fpara[6]/2
  #
  # 
  dis        = okadaDC3D(alpha,fx,fy,z,DEP,DIP,AL1,AL2,AW1,AW2,\
                         fpara[7],fpara[8],fpara[9]) 
  if strainonly:
      #
      # Strain tensor only
      # force sxy == syx
      # 
      xx,yy,zz,xy,xz,yz = dc3dstrain( dis[3]/1000.,(dis[4]+dis[6])/2000.,\
                                      (dis[5]+dis[9])/2000.,\
                                      dis[7]/1000.,(dis[8]+dis[10])/2000.,\
                                      dis[11]/1000,strike)
  else:
      xx,yy,zz,xy,xz,yz = dc3dstrain2stress(dis[3], dis[4], dis[5],\
                                            dis[6], dis[7], dis[8],\
                                            dis[9], dis[10],dis[11],\
                                            strike,young=young,pois=pois)
  #
  return xx,yy,zz,xy,xz,yz
#
###############################################################################
#
# @numba.jit(nopython=True, parallel=True, nogil=True)  
def dc3dstrain(uxx,uxy,uxz,uyy,uyz,uzz,\
                      strike):
    #
    strainmatrix = np.vstack((uxx,uyy,uzz,uxy,uxz,uyz)).T
    outstraintensor = stresstensorrot(strike,strainmatrix)
    sxx = outstraintensor[:,0]
    syy = outstraintensor[:,1]
    szz = outstraintensor[:,2]
    sxy = outstraintensor[:,3]
    sxz = outstraintensor[:,4]
    syz = outstraintensor[:,5]
    return sxx,syy,szz,sxy,sxz,syz
#
#
# @numba.jit(nopython=True, parallel=True, nogil=True)    
def dc3dstrain2stress(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,\
                      strike,young=1e6,pois=0.25):
    #
    # young, axial strain (elongations or shortenings)
    # It is fixed in default, based on a classic paper, as 1e6 bar, approximately 
    #    1e5 MPa
                      
    # Watts, A.B., 1983. The Strength of the Earth's Crust (No. LDGO-3449). 
    #                  LAMONT-DOHERTY GEOLOGICAL OBSERVATORY PALISADES NY.                 
    #                   
    #
    # shear modulus, sk
    # 
    # sk   = young/(2*(1.0+pois))
    # the first lame constant : lame1
    # http://scienceworld.wolfram.com/physics/LameConstants.html
    #
    # lame1 = (young * pois)/((1+pois)*(1-pois))
    # the second lame constant : lame2
    # lame2 = young / (2+2*pois)
    sk  = young/(1.0+pois)
    gk  = pois/(1.0-2.0*pois)
    # the trace of strain tensor
    vol = uxx + uyy + uzz 
    #
    # % caution! strain dimension is from x,y,z coordinate (should be /1000)
    #
    tsxx = sk * (gk * vol + uxx) * 0.001
    tsyy = sk * (gk * vol + uyy) * 0.001
    tszz = sk * (gk * vol + uzz) * 0.001
    tsxy = (young/(2.0*(1.0+pois))) * (uxy + uxy) * 0.001
    tsxz = (young/(2.0*(1.0+pois))) * (uxz + uxz) * 0.001
    tsyz = (young/(2.0*(1.0+pois))) * (uyz + uzy) * 0.001
    #
    stressmatrix = np.vstack((tsxx,tsyy,tszz,tsxy,tsxz,tsyz)).T
    #
    outstresstensor = stresstensorrot(strike,stressmatrix)
    #
    sxx = outstresstensor[:,0]
    syy = outstresstensor[:,1]
    szz = outstresstensor[:,2]
    sxy = outstresstensor[:,3]
    sxz = outstresstensor[:,4]
    syz = outstresstensor[:,5]
    return sxx,syy,szz,sxy,sxz,syz
############################################################################### 
def dc3dfromfaults(fpara,x,y,z,alpha=0.666666666667):
    #
    for i in range(fpara.shape[0]):
        ux,uy,uz = dc3dtodis(fpara[i,:],x,y,z,alpha=alpha)
        if i == 0:
            outx = np.copy(ux)
            outy = np.copy(uy)
            outz = np.copy(uz)
        else:
            outx = outx + ux
            outy = outy + uy
            outz = outz + uz
    return outx,outy,outz
    
#
# @numba.jit(nopython=True, parallel=True, nogil=True)    
def dc3dtodis(fpara,x,y,z,alpha=0.66666666666666667):
  #
  # Note that alpha used for okada3D is different with that for okada2D
  # 
  if not isinstance(x,np.ndarray):
     #
     # python version dc3d supports matrix operation, which can largely make 
     # computation faster. The input coordinates, (x,y,z) all should be 
     # in numpy.ndarray. If it was found not numpy.array, 
     # here we will force them to be ...
     # by Wanpeng Feng, @Ottawa, 2016-10-02
     #
     x  = np.array([x])
     y  = np.array([y])
     z  = np.array([z])
  #
  fpara      = np.array(fpara)
  fpara[3]   = (fpara[3] == 90.) * 89.9999 + (fpara[3] != 90. ) * fpara[3]
  fx,fy      = geo2okada3d(x,y,fpara)
  #
  # % alpha is media parameter. 
  diprad     = fpara[3]/180.*np.pi
  DEP        = fpara[4] + fpara[5]/2 * np.sin(diprad)
  strike     = fpara[2];
  DIP        = fpara[3]
  AW1        = fpara[5]/2 * -1
  AW2        = fpara[5]/2
  AL1        = fpara[6]/2 * -1
  AL2        = fpara[6]/2
  #
  dis        = okadaDC3D(alpha,fx,fy,z,DEP,DIP,AL1,AL2,AW1,AW2,\
                         fpara[7],fpara[8],fpara[9]) 
  #
  strkr      = np.complex(0,(90-strike)*np.pi/180)                
  disxy      = np.array(dis[0]*0,dtype='complex')
  disxy.real = np.copy(dis[0])
  disxy.imag = np.copy(dis[1])
  #
  disxy      = disxy * np.exp(strkr) 
  # dc3dstrain2stress(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,\
  #                    strike,young=800000,pois=0.25)
  return disxy.real, disxy.imag, dis[2]

# 
###############################################################################
#
def poisson2alp(nu,young=7.5*10e10,dens=2850):
    #
    # nu, poisson ratio
    #
    # dens is in k in m^3
    #
    # shear modulus
    #
    g = young/(2*(1+nu))
    m = 2*g*(1-nu)/(1-2*nu)
    vp = (m/dens)**(1./2)/1000
    vs = vp/((2*(1-nu)/(1-2*nu))**(1./2))
    alp = 1./((vp/vs)**2.-1)
    #
    return alp
#    
###############################################################################
def slipINV(fpara,inp,rake=[0,90],alpha=0.5,disp=False):
    '''
    To invert the data (e.g. InSAR, GPS and any other geodetic observations)
    for the slip with a linear least square algorithm. 
    Developed by Wanpeng Feng, @CCRS/NRCan, Canada, 2018-02-07
    '''
    mfpara = np.copy(fpara)
    G1,G2,rake = dc2dTOgreen(mfpara,inp,rake=rake,alpha=alpha,disp=disp)
    # output of rake through dc2dTOgreen should be in np.ndarray 
    #
    if sum(np.abs(rake[:,0]-rake[:,1])) > 0:
        G = np.hstack((G1,G2))
        flag = True
    else:
        G = G1
        flag = False
    #
    # slip the result 
    # res  the residual
    # rnk  Rank of matrix
    # s    singularity of G
    #
    slip,res,rnk,s = np.linalg.lstsq(G,inp[:,2],rcond=None)
    sim = np.dot(G,slip)
    #
    if flag:
       s1 = slip[0:mfpara.shape[0]]
       s2 = slip[mfpara.shape[0]:]
       str_slip = s1 * np.cos(rake[:,0]*np.pi/180.) + \
                  s2 * np.cos(rake[:,1]*np.pi/180.)
       dip_slip = s1 * np.sin(rake[:,0]*np.pi/180.) + \
                  s2 * np.sin(rake[:,1]*np.pi/180.)
    else:
       s1 = slip
       str_slip = s1 * np.cos(rake[:,0]*np.pi/180.)
       dip_slip = s1 * np.sin(rake[:,0]*np.pi/180.)
    #
    mfpara[:,7] = str_slip
    mfpara[:,8] = dip_slip
    return mfpara,res,rnk,s,sim
       
#
def dc2dTOgreen(mfpara,inp,rake=[0,90],alpha=0.5,disp=False):
    '''
    Calculate theoretical surface response of a unit slip on faults
    
    '''
    nf = mfpara.shape[0]
    #
    if rake is not None:
       #
       # Make sure rake is a list or an array
       #
       if isinstance(rake,(list,np.ndarray)):
          rake = np.array(rake)
          if len(rake.shape) == 1 and len(rake) == 2:
              rake = np.reshape(rake,[1,2])
       else:
          rake = np.reshape(np.array([rake,rake]),[1,2])
       #
       if rake.shape[0] < nf:
           rake = npm.repmat(rake[0:1,:],nf,1)
    else:
       rake = npm.repmat(np.reshape(np.array([0,90]),[1,2]),nf,1)
    #
    str_slip = np.cos(rake[:,0]*np.pi/180.) + np.cos(rake[:,1]*np.pi/180.)
    dip_slip = np.sin(rake[:,0]*np.pi/180.) + np.sin(rake[:,1]*np.pi/180.)
    #
    green_1 = np.zeros([inp.shape[0],nf])
    green_2 = np.copy(green_1)
    for i in range(nf):
        cfpara = mfpara[i:i+1,:]
        cfpara[0,7] = str_slip[i]
        cfpara[0,8] = 0
        sUe,sUn,sUv = dc2dfromfaults(cfpara,inp[:,0],inp[:,1],\
                                  alpha=alpha,disp=False)
        green_1[:,i] = sUe * inp[:,3] + sUn * inp[:,4] + sUv * inp[:,5]
        #
        cfpara[0,8] = dip_slip[i]
        cfpara[0,7] = 0
        dUe,dUn,dUv = dc2dfromfaults(cfpara,inp[:,0],inp[:,1],\
                                  alpha=alpha,disp=False)
        green_2[:,i] = dUe * inp[:,3] + dUn * inp[:,4] + dUv * inp[:,5]
    #
    return green_1,green_2,rake
       
#
# @numba.jit(nopython=True, parallel=True, nogil=True)    
def dc2dTOLOS(mfpara,inp,alpha=0.5,disp=False):
    '''
    To calculate the displacements along the light of sight for InSAR
    <inp>  an nx7 array
    output:
    dSIM   prediction of the data with a given model
    dRES   residuals of the data
    dSTD   standard deviations of residuals
    dMEAN  mean of residuals
    by Wanpeng Feng, @Ottawa, 2017-12-12
    
    Updated by FWP, @Ottawa, 2017-12-20
    dSIM will be returned as well.
    
    '''
    #
    Ue,Un,Uv = dc2dfromfaults(mfpara,inp[:,0],inp[:,1],alpha=alpha,disp=False)
    dSIM     = Ue * inp[:,3] + Un * inp[:,4] + Uv * inp[:,5]
    dRES     = inp[:,2] - dSIM
    dMEAN    = np.mean(dRES)
    dSTD     = np.std(dRES)
    #
    return dSIM,dRES,dMEAN,dSTD

#   
def dc2dfromfaults(mfpara,x,y,alpha=0.5,disp=False):
    #
    if len(mfpara) == 10:
        mfpara = np.reshape(mfpara,[1,10])
    #
    for i in range(mfpara.shape[0]):
        if disp:
          nofault = i+1
          print(' pokada: No: %-4d fault' % nofault)
        #
        xd,yd,vd = dc2dtodis(mfpara[i,:],x,y,alpha=alpha)
        #
        if i == 0:  
            outE = np.copy(xd)
            outN = np.copy(yd)
            outU = np.copy(vd)
        else:
            outE = outE + xd
            outN = outN + yd
            outU = outU + vd
            
    return outE,outN,outU

#          
def dc2dtodis(fpara,x,y,alpha=0.5):
  #
  fpara    = np.array(fpara)
  #
  # A sigularity for dip 90, to avoid a failure during calculation. the dip of 90
  # will be forced to 89.9999
  # by Wanpeng Feng, @Ottawa, 2016-10-03
  #
  fpara[3] = (fpara[3] == 90.) * 89.9999  + (fpara[3] != 90.) * fpara[3]
  # To avoid sigularity for points on the fault trace...
  fpara[4] = (fpara[4] == 0. ) * 0.00000001 + (fpara[4] != 0. ) * fpara[4]
  #
  fx,fy    = geo2okada2d(x,y,fpara)
  #
  # % alpha is media parameter. 
  SD     = np.sin(fpara[3]/180*np.pi)
  CD     = np.cos(fpara[3]/180*np.pi)
  DEP    = fpara[4] + fpara[5] * SD
  strike = fpara[2];
  #
  # print(SD,CD)
  dis        = okadaDC2D(alpha,fx,fy,DEP,0.,fpara[6],0.,fpara[5],\
                         SD,CD,fpara[7],fpara[8],fpara[9])    
  #
  strkr      = np.complex(0,(90-strike)*np.pi/180)                
  disxy      = np.array(dis[0],dtype='complex')
  disxy.real = np.copy(dis[0])
  disxy.imag = np.copy(dis[1])
  #
  disxy      = disxy *  np.exp(strkr) 
  
  return disxy.real, disxy.imag,dis[2]
  #
###############################################################################
def geo2okada3d(x,y,faultp):
  '''
  A python version of coordinate conversion was made by Wanpeng Feng based on my 
  previous PSOKINV inversion package. 
  This can be used to temporarily compute new coordinates with respect to a specific fault.
  This is a necesary step before running okada dislocation model for surface deformation calculation
  First done by Wanpeng Feng, @Ottawa, Candada, 2016-10-03
  #
  # #function [ox,oy] = geo2okada(x,y,faultp)
  '''
  xstart  = faultp[0]
  ystart  = faultp[1]
  dip     = faultp[3]
  strike  = faultp[2]
  width   = faultp[5]
  flen    = faultp[6]
  #
  strkr   = (90.-strike)*np.pi/180.
  strkr   = np.complex(0,strkr)
  raddip  = dip/180.*np.pi
  #############################################################################
  ## type 3, centre of fault 
  #fll     = np.complex(-flen/2,-width/2*np.cos(raddip))
  #flu     = np.complex(-flen/2, width/2*np.cos(raddip))
  #fru     = np.complex( flen/2, width/2*np.cos(raddip))
  #frl     = np.complex( flen/2,-width/2*np.cos(raddip))
  #
  # type 0, top-centre of fault 
  flla    = np.complex(-flen/2, -width*np.cos(raddip))
  flua    = np.complex(-flen/2,  0.)
  frua    = np.complex( flen/2,  0.)
  frla    = np.complex( flen/2, -width*np.cos(raddip))
  #
  fll     = np.complex(xstart,ystart) + flla*np.exp(strkr)
  flu     = np.complex(xstart,ystart) + flua*np.exp(strkr)
  frl     = np.complex(xstart,ystart) + frla*np.exp(strkr)
  fru     = np.complex(xstart,ystart) + frua*np.exp(strkr)
  #
  xstartd = (fll.real + flu.real + frl.real + fru.real) / 4
  ystartd = (fll.imag + flu.imag + frl.imag + fru.imag) / 4
  #############################################################################
  # 
  ccd     = np.complex(xstartd,ystartd) # + ll*np.exp(strkr);
  comxy   = np.array(x,dtype='complex')
  #
  comxy.real = np.copy(x)
  comxy.imag = np.copy(y)
  pokadaxy   = (comxy - ccd) * np.exp(-1.*strkr)
  #
  ox         = pokadaxy.real
  oy         = pokadaxy.imag
  # print(" Coverted COOR: %f %f " % (ox,oy))
  return ox,oy
###############################################################################  
# @numba.jit(parallel=True, nogil=True) 
def geo2okada2d(x,y,faultp):
  # #function [ox,oy] = geo2okada(x,y,faultp)
  xstart  = faultp[0]
  ystart  = faultp[1]
  dip     = faultp[3]
  strike  = faultp[2]
  width   = faultp[5]
  flen    = faultp[6]
  dip     = (dip==90)*89.9999 + (dip!=90.)*dip;
  ll      = np.complex(-1.*flen/2,-width*np.cos(dip/180.*np.pi)) 
  strkr   = (90.-strike)*np.pi/180.
  strkr   = np.complex(0,strkr)
  lld     = np.complex(xstart,ystart) + ll*np.exp(strkr);
  comxy   = np.array(x,dtype='complex')
  #
  comxy.real = np.copy(x)
  comxy.imag = np.copy(y)
  pokadaxy   = (comxy - lld) * np.exp(-1.*strkr)
  #
  ox         = pokadaxy.real
  oy         = pokadaxy.imag
  # print(" Coverted COOR: %f %f " % (ox,oy))
  return ox,oy
###############################################################################  
# @numba.jit(nopython=True,parallel=True, nogil=True)
def okadaDC2D(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,\
           SD,CD,DISL1,DISL2,DISL3):
  #
  # Transferend directly from okada85.f by wanpneg Feng, @Ottawa, 2016-10-01
  #            
  # function  OU=okada2d(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,...
  #                         SD,CD,DISL1,DISL2,DISL3)
  # %C*********************************************************  
  # %C*****                                               *****
  # %C*****    SURFACE DISPLACEMENT,STRAIN,TILT           *****
  # %C*****    DUE TO RECTANGULAR FAULT IN A HALF-SPACE   *****
  # %C*****       CODED BY Y.OKADA ... JAN 1985           ***** 
  # %C*****                                               *****
  # %C*********************************************************  
  # %C
  # %C***** INPUT
  # %C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)
  # %C*****   AL1,AL2 : FAULT LENGTH RANGE, a line vector#
  # %C*****   AW1,AW2 : FAULT WIDTH RANGE
  # %C*****   SD,CD   : SIN,COS OF DIP-ANGLE
  # %C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
  # %C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
  # %C
  # %C***** OUTPUT
  # %C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )
  # %C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
  # %C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )
  # % disp(DISL2);
  # % disp(DISL1);
  #         
  F0 = 0.;
  F1 = 1.0;
  #%C-----
  U = [[],[],[]]
  # %for i=1:9,U{i}=X.*0; end
  #%
  for i in range(3):
      U[i] = X * 0
      #
  P = Y * CD + DEP * SD #;   % Calculation P, size as X & Y
  Q = Y * SD - DEP * CD #;   % Calculation Q, size as X & Y

  for k in range(1,3):
    if k == 1:
        ET = P-AW1
    else:
        ET = P-AW2
    for j in range(1,3):
      if j == 1:
          XI = X - AL1
      else:
          XI = X - AL2
      jk = j+k
      if jk != 3: 
          SIGN = F1
      else:
          SIGN =-F1
      # print(type(XI),type(F0))
      XI = (XI!=F0) * XI + (XI==F0) * 0.0001
      ET = (ET!=F0) * ET + (ET==F0) * 0.0001
      #
      DU = MSRECTG_OKADA(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3);
      for i in range(3):
          U[i] = U[i] + SIGN * DU[i]
   
  #
  #
  return U
#
#%
##############################################################################
# @numba.jit(nopython=True,parallel=True, nogil=True)  
#function DU=MSRECTG_OKADA(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3)
def MSRECTG_OKADA(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3):
  # %C*********************************************************************
  # %C*****                                                           *****
  # %C*****  DUE TO RECTANGULAR FAULT IN A HALF-SPACE                 *****
  # %C*****                          CODED BY  Y.OKADA ... JAN 1985   *****  
  # %C*****                                                           *****
  # %C*********************************************************************
  # %C
  # %C***** INPUT
  # %C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)
  # %C*****   XI,ET,Q : FAULT COORDINATE
  # %C*****   SD,CD   : SIN,COS OF DIP-ANGLE
  # %C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
  # %C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
  # %C
  # %C***** OUTPUT
  # %C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )
  # %C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
  # %C*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )
  # %C
  F0  = 0.
  F1  = 1.
  F2  = 2.
  PI2 = 6.283185307179586
  #%
  XI2 = XI ** 2
  ET2 = ET ** 2
  Q2  = Q  ** 2
  R2  = XI2 + ET2 + Q2
  R   = np.sqrt(R2)
  # R3  = R * R2
  D   = ET * SD - Q * CD
  Y   = ET * CD + Q * SD
  RET = R  + ET
  #
  if isinstance(RET,(np.ndarray)):
      RET[RET<F0] = F0 
  else:
      if RET < F0:
          RET = F0;
  #RET[RET< F0] = F0
  RD  = R+D
  #%RRD=F1./(R.*RD);
  #%
  TT  = (Q != F0) * np.arctan( XI*ET/(Q*R)) + (Q==F0) * F0
  #%
  RE  = (RET != F0) * F1/RET   + (RET==F0) * F0
  DLE = (RET != F0) * np.log(RET) - (RET==F0) * np.log(R-ET)
  RRX = F1 / (R * (R+XI))
  RRE = RE / R
  if CD == F0: 
     RD2 = RD*RD
     A1  = -ALP /F2 * XI * Q / RD2 
     A3  =  ALP /F2 * ( ET /RD + Y * Q / RD2 - DLE ) 
     A4  = -ALP * Q / RD;
     A5  = -ALP * XI * SD / RD 
  else:
     #%C==============================
     #%C=====   INCLINED FAULT   =====
     #%C==============================
     TD = SD / CD
     X  = np.sqrt(XI2 + Q2)
     # %A5=(XI==F0).*XI.*F0+(XI~=F0).*XI.*ALP.*F2./CD.*atan((ET.*(X+Q.*CD)+X.*(R+X).*SD)./(XI.*(R+X).*CD));
     A5 = (XI==F0) * F0 + (XI != F0) * ALP * F2 / CD * np.arctan((ET*(X+Q*CD)+X*(R+X)*SD)/((R+X)*XI*CD)) 
     A4 = ALP / CD * (np.log(RD) - SD * DLE)
     A3 = ALP * (Y/RD/CD - DLE) + TD * A4
     A1 = -ALP/CD*XI/RD - TD*A5
  #%
  A2 = -ALP * DLE-A3
  U1 =  XI * 0
  U2 =  XI * 0
  U3 =  XI * 0
  # %C======================================
  # %C=====  STRIKE-SLIP CONTRIBUTION  =====
  # %C======================================
  UN  = (DISL1 != F0) * DISL1 / PI2 + (DISL1 == F0) * 0   
  REQ = RRE * Q
  U1  = U1 - UN * ( REQ * XI + TT          + A1 * SD )
  U2  = U2 - UN * ( REQ * Y  + Q * CD * RE + A2 * SD )
  U3  = U3 - UN * ( REQ * D  + Q * SD * RE + A4 * SD )
  # %C===================================
  # %C=====  DIP-SLIP CONTRIBUTION  =====
  # %C===================================
  UN   = (DISL2 != F0) * DISL2 / PI2 + (DISL2 == F0) * 0
  SDCD = SD * CD
  U1   = U1 - UN * ( Q / R - A3 * SDCD );
  U2   = U2 - UN * ( Y * Q * RRX + CD * TT - A1 * SDCD )
  U3   = U3 - UN * ( D * Q * RRX + SD * TT - A5 * SDCD )
  # %C=========================================================================
  # %C===================  TENSILE-FAULT CONTRIBUTION  ========================
  # %C=========================================================================
  UN   = (DISL3 != F0) * DISL3 / PI2 + (DISL3==F0) *0
  SDSD = SD * SD
  U1   = U1 + UN * ( Q2 *RRE - A3 * SDSD )
  U2   = U2 + UN * (-D * Q * RRX - SD * (XI * Q * RRE - TT) - A1 * SDSD )
  U3   = U3 + UN * ( Y * Q * RRX + CD * (XI * Q * RRE - TT) - A5 * SDSD )
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DU    = [0,0,0]  
  DU[0] = U1 #%*FACTOR;
  DU[1] = U2 #%*FACTOR;
  DU[2] = U3 #%*FACTOR; 
  
  return DU
 
