import numpy as np
import matplotlib.pyplot as plt


def JSWP(HS,TP,GAMMA,FREQ):
 nfreq=len(FREQ)
# SPECTRE DE JONSWAP
 DELTA=0.0624/(0.23+0.033*GAMMA-0.185/(1.9+GAMMA))*HS**2
 FP=1.0/TP
 DELTA=DELTA*FP**4
 SIGMAA=0.07000
 SIGMAB=0.09000
 TMIN=TP/3.0
 TMAX=2.0*TP
 FMIN=FREQ[0]
 FMAX=FREQ[-1]

 JONSWAP=[]
 for i in range(nfreq):
       FF=FREQ[i]
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
# here we could do something better (case when df is not the same for any f)
 df  =float(FMAX-FMIN)/float(nfreq)
 AMP    = np.sqrt(2.0*df*JONSWAP) 
 return AMP



