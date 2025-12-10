import numpy as np
import matplotlib.pyplot as plt



def readRad(filename):
  f1=open(filename,'r')
  text1=f1.readlines()
  nfreq=(len(text1)-1)//7
  CM      =np.zeros((nfreq,1+36))
  MatOmega=np.zeros((nfreq,3,3))
  # lecture full file
  for i in range(nfreq):
     omeg=[float(x) for x in text1[1+7*i].split()]
     CM[i,0]=omeg[0]
     for il in range(6):
      lfloats=[float(x) for x in text1[1+7*i+il+1].split()]
      CM[i,1+il*6:1+il*6+6]=lfloats
  for i in range(nfreq):
   MatOmega[i,0,0]=CM[i,1+0*6+0]
   MatOmega[i,1,1]=CM[i,1+2*6+2]
   MatOmega[i,2,2]=CM[i,1+4*6+4]
   MatOmega[i,0,2]=CM[i,1+0*6+4]
   MatOmega[i,2,0]=CM[i,1+4*6+0]
  # keep only DoF of interest
  return MatOmega


def cutDDL(mat66):
  mat33=np.zeros((3,3))
  mat33[0,0]=mat66[0,0]
  mat33[1,1]=mat66[2,2]
  mat33[2,2]=mat66[4,4]
  mat33[0,2]=mat66[0,4]
  mat33[2,0]=mat66[4,0]
  return mat33



def cutRAO(rao6):
  nfreq=len(rao6[:,0])
  rao3=np.zeros((nfreq,1+6))
  for ifreq in range(nfreq):
    rao3[ifreq,0]=rao6[ifreq,0]
    rao3[ifreq,1]=rao6[ifreq,1]
    rao3[ifreq,2]=rao6[ifreq,3]
    rao3[ifreq,3]=rao6[ifreq,5]
    rao3[ifreq,4]=rao6[ifreq,7]
    rao3[ifreq,5]=rao6[ifreq,9]
    rao3[ifreq,6]=rao6[ifreq,11]
  return rao3




def force_wind(Ws,D,Cd,HH,iddl):
  F=0.
  coeff=0.
  if iddl==0:
    coeff=1.
  if iddl==1:
    coeff=0.
  if iddl==2:
    coeff=HH
  F=0.5*Cd*1.2*np.pi*D**2/4. * Ws**2 * coeff
  return F

