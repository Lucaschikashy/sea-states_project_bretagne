import numpy as np
import matplotlib.pyplot as plt


###############################################################################################################################
## Sea State
# Hauteur significative
HS       = 10.  #m
# Periode Pic
TP       = 16.1  #s
# Facteur de pic
GAMMA = 1.
# Duree de l'etat de mer
Duration = 10800.0 #s 10800s <=> 3h
# Nombre de frequences dans le spectre
NE       = 100 
###############################################################################################################################
















###############################################################################################################################


# SPECTRE DE JONSWAP
DELTA=0.0624/(0.23+0.033*GAMMA-0.185/(1.9+GAMMA))*HS**2
FP=1.0/TP
DELTA=DELTA*FP**4
SIGMAA=0.07000
SIGMAB=0.09000
TMIN=TP/3.0
TMAX=2.0*TP
FMIN=1.0/TMAX
FMAX=1.0/TMIN
df  =float(FMAX-FMIN)/float(NE)


dt= 0.1
N = int(Duration/dt)


FREQJ=[]
JONSWAP=[]
for i in range(NE):
       FREQJ.append(FMIN+(i-1)*(FMAX-FMIN)/NE)
       FF=FMIN+(i-1)*(FMAX-FMIN)/NE
       if (FF<FP) :
         SIG=SIGMAA
       else:
         SIG=SIGMAB

       ARG1=0.50*((FF-FP)/(SIG*FP))**2
       if (ARG1<99.0) :
         ARG1=GAMMA**np.exp(-ARG1)
       else:
         ARG1=1.0

       ARG2=1.250*(FP/FF)**4
       if (ARG2<99.0) :
         ARG2=np.exp(-ARG2)
       else:
         ARG2=0.
 
       ARG3=DELTA/FF**5
       JONSWAP.append(ARG1*ARG2*ARG3)

JONSWAP=np.array(JONSWAP)
AMP    = np.sqrt(2.0*df*JONSWAP) 
PHAS   = np.random.rand(NE) * 2.0*np.pi

FREQJ  =np.array(FREQJ)
OMEGJ  = 2.0*np.pi * FREQJ
OMEGJ  = OMEGJ + 0.5*PHAS*df


Time = np.zeros(N)
Eta  = np.zeros(N)
for i in range(N) :
  Time[i] = dt*float(i)
  CPHASJ  = np.cos(OMEGJ*Time[i] + PHAS) 
  Eta[i]  = np.dot(AMP,CPHASJ) 
###############################################################################################################################



plt.plot(Time , Eta,"-b",label="Surface Libre")
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Surface Libre (m)')
plt.legend()  
plt.show()



