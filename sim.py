import matplotlib.pyplot as plt
import numpy as np
from math import *

#convenient functions
sign=lambda v:abs(v)/v if v!=0 else 0
def npmap(l, f):
    return np.array(list(map(f, list(l))))
def Butcher_Tableaux(a,c,b,N=20,Err=.00001):
    k=[lambda f,t,dt,u:0 for i in range(len(a))]
    csum=[lambda f,t,dt,u:0 for i in range(len(a))]
    for i in range(len(a)):
        print(i)
        if sum(a[i][i+1:])!=0 :
            raise Exception("tableaux a has values above the diagonal")
        print([[j,a[i][j]] for j in range(i)])
        csum[i]=lambda f,t,dt,u:sum([a[i][j]*k[j](f,t,dt,u) for j in range(i)])
        if a[i][i]!=0:#check if implicit
            k[i]=lambda f,t,dt,u: iterate(lambda un:f(t+c[i]*dt,u[-1]
                +dt*(csum[i](f,t,dt,u)+a[i][i]*un)),f(t,u[-1]))
        else:#it's explicit :)
            k[i]=lambda f,t,dt,u: f(t+c[i]*dt,u[-1]+dt*csum[i](f,t,dt,u))
        return lambda f,t,dt,u:(u[-1]
            +dt*sum([b[j]*k[j](f,t,dt,u) for j in range(len(b))]))

#given a function itertate f such that f(f(f(...))) converges find that 
def iterate(f,guess,iter=20,err=0.00001):
    i=0
    fguess=f(guess)
    while(i<iter and np.linalg.norm(fguess-guess)>err):
        guess=fguess
        fguess=f(guess)
        i+=1
    return guess
def implicit(iterable):
    return lambda f,t,dt,u: iterate(iterable(f,t,dt,u),u[-1])
#Finite Difference Scheme(FDS) integrator
def integrate(u0, scheme, f, T, dt):
    u = u0
    print(u)
    N=ceil(T/dt)
    for n in range(N):
        u+=[scheme(f,dt*n,dt,u)]
        if u[-1][0]<0:
            return u
    return u

#initial conditions
T=200;dt=0.05;
U0=np.array([0,0])#u[0]=height u[1]=velocity

#physical constants
def DrafCoeff(h,v):#TODO consider Angle
    return -sign(v)*0.55
def AirDensity(h):
    return 1.2*.99988**h
def Thrust(t):
    return 3000 if t<4 else 0
g=-10
Area=(pi/4)*(6/39)**2#TODO consider angle
M=25#TODO consider change in mass

#acceleration=-(A*rho*Cd*v^2+thrust)/m+g
Accel=lambda t,h,v:(Area*AirDensity(h)*DrafCoeff(h,v)*v**2+Thrust(t))/M+g

#FDS bs
def F(t,u):#the derivative of the state space
    return np.array([u[1],Accel(t,u[0],u[1])])
def FE(f,t,dt,u):#forward euler FDS
    return u[-1]+dt*f(t,u[-1])
def BE(f,t,dt,u):#Backward Euler
    return iterate(lambda un:u[-1]+dt*f(t+dt,un),u[-1])
def Trap(f,t,dt,u):#trapezoidal FDS
    return iterate(lambda un:u[-1]+dt*(f(t+dt,un)+f(t,u[-1]))/2,u[-1])
RK4=Butcher_Tableaux([[0   , 0,0,0],\
                      [ 1/3, 0,0,0],\
                      [-1/3, 1,0,0],\
                      [1   ,-1,1,0]],[0,1/3,2/3,1],[1/8,3/8,3/8,1/8])
RK1=Butcher_Tableaux([[0]],[0],[1])

#run code
U=np.array(integrate([U0,U0],FE,F,T,dt))
UBE=np.array(integrate([U0,U0],RK1,F,T,dt))
UTrap=np.array(integrate([U0,U0],Trap,F,T,dt))

# plotting
plt.plot(dt*np.arange(len(UTrap)),UTrap[:,0],label="Trapy")
plt.plot(dt*np.arange(len(UTrap)),UTrap[:,1],label="Trapv")
plt.plot(dt*np.arange(len(UBE)),UBE[:,0],label="BEy")
plt.plot(dt*np.arange(len(UBE)),UBE[:,1],label="BEv")
plt.plot(dt*np.arange(len(U)),U[:,0],label="FEy")
plt.plot(dt*np.arange(len(U)),U[:,1],label="FEv")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("velocity(m/s)/Height(m)")
plt.show()
