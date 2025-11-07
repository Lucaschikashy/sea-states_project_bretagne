import numpy as np
import matplotlib.pyplot as plt
from ocean_wave_tracing import Wave_tracing
import cmocean
from pathlib import Path

period= 5.0 # wave period (s)
ang   = 105 # wave direction (deg)
nr    = 30  # number of wave rays
simtime = 3600*1.5 # propagate 1.5 hour
nt    = 1200 # number of time steps
fname = Path(__file__).resolve().parent / 'GMRTv4_4_0_20251107topo.asc' # bathymetry file

# bathymetry
fid = open(fname)
ncols = int(fid.readline().split()[1])
nrows = int(fid.readline().split()[1])
lonll = float(fid.readline().split()[1])
latll = float(fid.readline().split()[1])
dd = float(fid.readline().split()[1])
dx = dd*6378.e3*np.pi/180*np.cos(np.pi/180*latll)
dy = dd*6378.e3*np.pi/180
fid.close()
bathy = -np.flipud(np.loadtxt(fname,skiprows=6))

# ray tracing
ang = 90-ang # convert
wt = Wave_tracing(U=np.zeros((nrows,ncols)),
                  V=np.zeros((nrows,ncols)),
                  nx=ncols, ny=nrows, nt=nt, T=simtime,
                  dx=dx,dy=dy,nb_wave_rays=nr,
                  domain_X0=0, domain_XN=dx*(ncols-1),
                  domain_Y0=0, domain_YN=dy*(nrows-1),
                  d=bathy)
if np.mod(ang,360)<45 or np.mod(ang,360)>=315:
    iws = 'left'
elif np.mod(ang,360)<135 and np.mod(ang,360)>=45:
    iws = 'bottom'
elif np.mod(ang,360)<225 and np.mod(ang,360)>=135:
    iws = 'right'
else:
    iws = 'top'
wt.set_initial_condition(wave_period=period,
                         theta0=+np.pi/180*ang,
                         incoming_wave_side=iws)
wt.solve()

# plotting wave rays
fig, ax = plt.subplots()
ax.set_aspect('equal')
cmap = cmocean.cm.deep.reversed()
cmap.set_bad(color='gray')
pc=ax.pcolormesh(wt.x,wt.y,-wt.d,shading='auto',cmap=cmap)
fig.colorbar(pc)
for ray_id in range(wt.nb_wave_rays):
    ax.plot(wt.ray_x[ray_id,:],wt.ray_y[ray_id,:],'-k')
ax.set_xlim((0,dx*(ncols-1)))
ax.set_ylim((0,dy*(nrows-1)))
pc.set_clim((-np.max(wt.d),0))
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Wave rays')
plt.show()

# wave energy flux density estimation
from scipy.interpolate import RegularGridInterpolator
xx,yy,hm=wt.ray_density(x_increment=20,y_increment=20)
interp = RegularGridInterpolator(points=(xx[0,:],yy[:,0]),values=hm.transpose(),bounds_error=False,fill_value=None)
xg,yg=np.meshgrid(wt.x,wt.y,indexing='ij')
e = interp((xg,yg))
e = e.transpose()
e[np.isnan(wt.d)]=np.nan
if np.mod(ang,360)<45 or np.mod(ang,360)>=315:
    escale = e[round(nrows/2),0]
elif np.mod(ang,360)<135 and np.mod(ang,360)>=45:
    escale = e[0,round(ncols/2)]
elif np.mod(ang,360)<225 and np.mod(ang,360)>=135:
    escale = e[round(nrows/2),-1]
else:
    escale = e[-1,round(ncols/2)]
fig, ax = plt.subplots()
ax.set_aspect('equal')
cmap = plt.cm.get_cmap("jet")
cmap.set_bad(color='gray')
pc = ax.pcolormesh(wt.x,wt.y,e/escale,shading='auto',cmap=cmap)
fig.colorbar(pc)
ax.set_xlim((0,dx*(ncols-1)))
ax.set_ylim((0,dy*(nrows-1)))
pc.set_clim((0,np.max(e[~np.isnan(e)]/escale)))
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Relative wave energy flux density')
plt.show()
