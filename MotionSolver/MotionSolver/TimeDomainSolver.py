import numpy as np
import matplotlib.pyplot as plt
import spectrum as SS
import floating_reader as FL
from scipy import signal 

###############################################################################################################################
# Significant Wave Height
Hs       = 0.0   #m
# Peak Period 
Tp       = 14.    #s
# Peak enhancement factor 
gamma    =  2.
# Number of wave components in the spectrum
nfreq    = 200
# Sea State duration
Duration = 1000.0 #s

# Wind speed 
Ws=0. #12.
# Turbine Diameter
Diam=178.
# Turbine height above sea level
HubHeight=120.
# Drag coefficient of the turbine
Cd=0.75

# Initial Surge
Surge0=0.
# Initial Heave
Heave0=0.
# Initial Pitch
Pitch0=10.

# Time step
dt= 0.1

###############################################################################################################################

# FLOATER
floater_folder='barge_data'

# Reading inputs files
MIfull =np.loadtxt('./'+floater_folder+'/MassInertia.dat',skiprows=2)
KMfull =np.loadtxt('./'+floater_folder+'/MooringStiffness.dat',skiprows=2)
KHfull =np.loadtxt('./'+floater_folder+'/HydrostaticStiffness.dat',skiprows=2)
BQfull =np.loadtxt('./'+floater_folder+'/QuadraticDamping.dat',skiprows=2)
RAOfull=np.loadtxt('./'+floater_folder+'/Fe.dat',skiprows=4)
CMfull =FL.readRad('./'+floater_folder+'/CM.dat')
CAfull =FL.readRad('./'+floater_folder+'/CA.dat')

# Set at python 3DoF format
MI =FL.cutDDL(MIfull)
KM =FL.cutDDL(KMfull)
KH =FL.cutDDL(KHfull)
BQ =FL.cutDDL(BQfull)
MA =np.zeros((3,3))
RD =np.zeros((3,3))
RAO=FL.cutRAO(RAOfull)


# Freq definition
FREQ=np.zeros(nfreq)
fmin=(1.0/Tp)*0.5
fmax=(1.0/Tp)*3.0
for i in range(nfreq):
  FREQ[i]=fmin +float(i)/float(nfreq-1)*(fmax-fmin)

# Outputs definition
N = int(Duration/dt)
Time = np.zeros(N)
Eta  = np.zeros(N)
X    = np.zeros((N,3))
XP   = np.zeros((N,3))
XPP  = np.zeros((N,3))
Faero=np.zeros((N,3))

# Initial position
X[0,0]=Surge0
X[0,1]=Heave0
X[0,2]=Pitch0*np.pi/180.

# Free Surface Amplitudes/Phases
AMP   =SS.JSWP(Hs,Tp,gamma,FREQ)
PHASE =np.random.rand(nfreq) * 2.0*np.pi

# Forces in the frequency domain: 
F=np.zeros((nfreq,3))
P=np.zeros((nfreq,3))
for ifreq in range(nfreq):
  for iddl in range(3):
     F[ifreq,iddl]=AMP[ifreq] *np.interp(2.0*np.pi*FREQ[ifreq],RAO[:,0],RAO[:,iddl+1])
     P[ifreq,iddl]=np.pi/180.0*np.interp(2.0*np.pi*FREQ[ifreq],RAO[:,0],RAO[:,iddl+1+3])

# Added mass & radiation damping
for ifreq in range(len(RAO[:,0]-2)):
  if (2.0*np.pi/Tp)<RAO[ifreq,0] and (2.0*np.pi/Tp)>RAO[ifreq-1,0]:
    iref=ifreq
coeff=1.0
if Hs<1e-6 :
  coeff=0.
  iref=0
MA=CMfull[iref,:,:]
RD=CAfull[iref,:,:]*coeff

###############################################   
# Time domain simulation
###############################################   
Force=np.zeros(3)
for i in range(N):
  Time[i] = dt*float(i)
  ramp = min(Time[i]/120.0,1.0)
  CPHASJ  = np.cos(-2.0*np.pi*FREQ*Time[i] + PHASE) 
  Eta[i]  = np.dot(AMP,CPHASJ)*ramp 
# Euler Scheme
  if i>0:
    X[i,:]    = X[i-1,:]+dt*XP[i-1,:]
    aXPXP       =abs(XP[i-1,:]) * XP[i-1,:]
    for iddl in range(3):
      Faero[i,iddl] = FL.force_wind(Ws,Diam,Cd,HubHeight,iddl)
      CPHASJ       = np.cos(-2.0*np.pi*FREQ*Time[i] + PHASE + P[:,iddl]) 
#                           Hydro Loads              +  Aero Loads
      Force[iddl]  =     (np.dot(F[:,iddl],CPHASJ)    + Faero[i,iddl]              )*ramp
    XPP[i,:]= np.dot(   np.linalg.inv(MI+MA), Force  -np.dot(RD,XP[i-1,:]) -np.dot(BQ,aXPXP) -np.dot(KM+KH,X[i,:])   )
    XP[i,:]=XP[i-1,:]+dt*XPP[i,:]
###############################################
###############################################   


###############################################
# Outputs in the time domain
###############################################   
plt.plot(Time , Eta,"-b",label="Surface Libre")
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Surface Libre (m)')
plt.legend()  
plt.show()

plt.subplot(311)
plt.plot(Time , X[:,0],"-b",label="Surge")
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Motion (m)')
plt.legend()  

plt.subplot(312)
plt.plot(Time , X[:,1],"-b",label="Heave")
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Motion (m)')
plt.legend()  

plt.subplot(313)
plt.plot(Time , X[:,2]*180.0/np.pi,"-b",label="Pitch")
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Motion (deg)')
plt.legend()  
plt.show()


###############################################
# Outputs in the frequency domain
###############################################   
fs1=1.0/dt
f1, Sp1 = signal.welch(Eta,fs1,window='hann',scaling='density',nperseg=int(N/10),noverlap=int(N/20),return_onesided=True)
plt.plot(f1,np.log10(Sp1), "-b",label='PSD Eta')
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlim([0.0,0.5])
plt.legend()
plt.show()


plt.subplot(311)
fs1=1.0/dt
f1, Sp1 = signal.welch(X[:,0],fs1,window='hann',scaling='density',nperseg=int(N/10),noverlap=int(N/20),return_onesided=True)
plt.plot(f1,np.log10(Sp1), "-b",label='PSD Surge')
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlim([0.0,0.5])
plt.legend()

plt.subplot(312)
fs1=1.0/dt
f1, Sp1 = signal.welch(X[:,1],fs1,window='hann',scaling='density',nperseg=int(N/10),noverlap=int(N/20),return_onesided=True)
plt.plot(f1,np.log10(Sp1), "-b",label='PSD Heave')
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlim([0.0,0.5])
plt.legend()

plt.subplot(313)
fs1=1.0/dt
f1, Sp1 = signal.welch(X[:,2],fs1,window='hann',scaling='density',nperseg=int(N/10),noverlap=int(N/20),return_onesided=True)
plt.plot(f1,np.log10(Sp1), "-b",label='PSD Pitch')
#plt.grid(b=True, which='major', color='k', linestyle='dotted',linewidth=1)
plt.xlim([0.0,0.5])
plt.legend()

plt.show()

###############################################
###############################################


