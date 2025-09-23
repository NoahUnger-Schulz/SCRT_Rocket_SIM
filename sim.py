import matplotlib.pyplot as plt
import numpy as np
from math import *

#convenient functions
sign=lambda v:abs(v)/v if v!=0 else 0
def npmap(l, f):
    return np.array(list(map(f, list(l))))

#Finite Difference Scheme(FDS) integrator
def integrate(u0, scheme, f, T, dt):
    u = u0
    N=ceil(T/dt)
    for n in range(N):
        u+=[scheme(f,dt*n,dt,u)]
        if u[-1][0]<0:
            return u
    return u

#initial conditions
T=200;dt=0.0005;
U0=np.array([0,0])#u[0]=height u[1]=velocity

#physical constants
DrafCoeff= lambda h,v: -sign(v)*0.55
AirDensity= lambda h:1.2*.99988**h
g=-10
Area=(pi/4)*(6/39)**2
Thrust=lambda t:3000 if t<4 else 0
M=20
Accel=lambda t,h,v:(Area*AirDensity(h)*DrafCoeff(h,v)*v**2+Thrust(t))/M+g

#FDS bs
F=lambda t,u: np.array([u[1],Accel(t,u[0],u[1])])
FE=lambda f,t,dt,u: u[-1]+dt*f(t,u[-1])
#CentralExp=lambda f,t,dt,u: 2*u[-1]-u[-2]-f(t,u[-1])*dt**2

#run code
U=np.array(integrate([U0,U0],FE,F,T,dt))

#plotting
plt.plot(dt*np.arange(len(U)),U[:,0])
plt.plot(dt*np.arange(len(U)),U[:,1])
plt.xlabel("time (s)")
plt.ylabel("velocity(m/s)/Height(m)")
plt.show()
