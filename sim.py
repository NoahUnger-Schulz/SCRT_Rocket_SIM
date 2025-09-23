import matplotlib.pyplot as plt
import numpy as np
from math import *

def npmap(l, f):
    return np.array(list(map(f, list(l))))


def integrate(u0, scheme, f, T, dt):
    u = u0
    N=ceil(T/dt)
    for n in range(N):
        u+=[scheme(f,dt*n,dt,u)]
        if u[-1][0]<0:
            return u
    return u

U0=np.array([0,0])
sign=lambda v:abs(v)/v if v!=0 else 0
Cd= lambda h,v: -sign(v)*0.55
rho= lambda h:.99988**h
g=-10
Thrust=lambda t:3000 if t<4 else 0
M=20

T=200;dt=0.0005;

F=lambda t,u: np.array([u[1],(((pi/4)*(6/39)**2)*(1.2)*rho(u[0])*Cd(u[0],u[1])*u[1]**2+Thrust(t))/M+g])

FE=lambda f,t,dt,u: u[-1]+dt*f(t,u[-1])

CentralExp=lambda f,t,dt,u: 2*u[-1]-u[-2]-f(t,u[-1])*dt**2


U=np.array(integrate([U0,U0],FE,F,T,dt))
print(np.array(U)[:,0])
plt.plot(dt*np.arange(len(U)),U[:,0])
plt.plot(dt*np.arange(len(U)),U[:,1])
plt.xlabel("time (s)")
plt.ylabel("velocity(m/s)/Height(m)")
plt.show()
