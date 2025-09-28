import matplotlib.pyplot as plt
import time
import numpy as np
from math import *

#convenient functions
sign=lambda v:abs(v)/v if v!=0 else 0
def npmap(l, f):
    return np.array(list(map(f, list(l))))
def Butcher_Tableaux(a,c,b,N=20,Err=.00001):
    k=[0 for i in range(len(a))]
    csum=[0 for i in range(len(a))]
    print(len(a),"len(a)")
    for i in range(len(a)):
        print(i)
        if sum(a[i][i+1:])!=0 :
            raise Exception("tableaux a has values above the diagonal")
        print([[j,a[i][j]] for j in range(i)])
        csum[i]=lambda f,t,dt,u,asum=csum[i-1]:sum([a[i][j]*k[j](f,t,dt,u,asum) for j in range(i)])
        if a[i][i]!=0:#check if implicit
            k[i]=lambda f,t,dt,u,asum=csum[i]: iterate(lambda un:f(t+c[i]*dt,u[-1]
                +dt*(asum(f,t+c[i]*dt,dt,u)+a[i][i]*un)),f(t,u[-1]))
        else:#it's explicit :)
            k[i]=lambda f,t,dt,u,asum=csum[i]: f(t+c[i]*dt,u[-1]+dt*asum(f,t,dt,u))
    return lambda f,t,dt,u:(u[-1]
            +dt*sum([b[j]*k[j](f,t,dt,u) for j in range(len(b))]))

def Butcher_Table(a,b,c):
    def k (f,t,dt,u,i):
        if i>=0:
            asum=np.array([0.,0.])
            for j in range(i):
                print(a[i][j],"aij")
                asum+=a[i][j]*k(f,t,dt,u,j)
            print(asum,"asum")
            if a[i][i]!=0:
                return iterate(lambda un:f(t+c[i]*dt,u[-1]
                    +dt*(asum+a[i][i]*un)),f(t+c[i]*dt,u[-1]))
            else:
                return f(t+c[i],u[-1]+dt*asum)
        else:
            return f(t+c[i]*dt,u[-1])
    return lambda f,t,dt,u:(u[-1]
            +dt*sum([b[j]*k(f,t,dt,u,j) for j in range(len(b))]))

def Butcher_Table2(a,b,c):
    def calc(f,t,dt,u):
        k=[]
        for i in range(len(c)):
            asum=np.array([0.,0.])
            for j in range(i):
                asum+=a[i][j]*k[j]
            k+=[f(t+c[i]*dt,u[-1]+asum*dt)]
            if a[i][i]!=0:
                k[-1]=iterate(lambda un:f(t+c[i]*dt,u[-1]
                        +dt*(asum+a[i][i]*un)),k[-1])
        bsum=np.array([0.,0.])
        for j in range(len(b)):
            bsum+=b[j]*k[j]
        return u[-1]+dt*bsum
    return calc

#given a function itertate f such that f(f(f(...))) converges find that 
def iterate(f,guess,iter=20,err=0.00001):
    i=0
    fguess=f(guess)
    while(i<iter and np.linalg.norm(fguess-guess)>err):
        guess=fguess
        fguess=f(guess)
        i+=1
    if i>1:
        print(i,"iters")
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

def double_step(scheme):
    return lambda f,t,dt,u:scheme(f,t+dt/2,dt/2,u+[scheme(f,t,dt/2,u)])

def integrate_adaptive(u0, scheme, f, T, dt,batch=1,scheme2=None,tol=None):
    if tol==None:
        tol=dt/1000
    if scheme2==None:
        scheme2=double_step(scheme)
    u = [u0[0] for i in range(batch+1)]
    t = [0 for i in range(batch+1)]
    print(u)
    while t[-1]<T:
        err=np.linalg.norm(scheme(f,t[-1],dt,u)-scheme2(f,t[-1],dt,u))
        #err/=np.linalg.norm(scheme(f,t[-1],dt,u))
        if err>tol:
            #print(err,dt,t[-1])
            u=u[:-batch]
            t=t[:-batch]
            dt/=4.0
        else:
            dt*=1.1
        for i in range(batch):
            t+=[t[-1]+dt]
            u+=[scheme(f,t[-1],dt,u)]
            if u[-1][0]<0 or t[-1]>T:
                return np.array(u),np.array(t)
    return np.array(u),np.array(t)





#initial conditions
T=200;dt=.05;
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
def Heun(f,t,dt,u):#forward euler FDS
    return u[-1]+dt*(f(t,u[-1])+f(t+dt,u[-1]+dt*f(t,u[-1])))/2
def FE(f,t,dt,u):#forward euler FDS
    return u[-1]+dt*f(t,u[-1])
def BE(f,t,dt,u):#Backward Euler
    return iterate(lambda un:u[-1]+dt*f(t+dt,un),u[-1])
def Trap(f,t,dt,u):#trapezoidal FDS
    return iterate(lambda un:u[-1]+dt*(f(t+dt,un)+f(t,u[-1]))/2,u[-1])

RK4=Butcher_Table2([[0   , 0,0,0],\
                   [ 1/3, 0,0,0],\
                   [-1/3, 1,0,0],\
                   [1   ,-1,1,0]],[0,1/3,2/3,1],[1/8,3/8,3/8,1/8])
RK1=Butcher_Table2([[1]],[1],[1])
RKtrap=Butcher_Table2([[0  ,0  ],\
                      [1/2,1/2]],[0,1],[0.5,0.5])
RKheun=Butcher_Table2([[0  ,0  ],\
                      [1,0]],[0,1],[0.5,0.5])
#run code
start=[]
start+=[(time.perf_counter())]
U,Time=integrate_adaptive([U0,U0],Heun,F,T,dt/10000)
start+=[(time.perf_counter())]
print(start[-1]-start[-2])
UTrap=np.array(integrate([U0,U0],FE,F,T,dt/100))
start+=[(time.perf_counter())]
print(start[-1]-start[-2])
UBE=np.array(integrate([U0,U0],FE,F,T,dt/1000))
start+=[(time.perf_counter())]
print(start[-1]-start[-2])

# plotting
plt.plot(dt/100*np.arange(len(UTrap)),UTrap[:,0],label="dt/100")
plt.plot(dt/100*np.arange(len(UTrap)),UTrap[:,1],label="dt/100")
plt.plot(dt/1000*np.arange(len(UBE)),UBE[:,0],label="dt/1000")
plt.plot(dt/1000*np.arange(len(UBE)),UBE[:,1],label="dt/1000")
plt.plot(Time,U[:,0],label="adaptive")
plt.plot(Time,U[:,1],label="adaptive")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("velocity(m/s)/Height(m)")
plt.show()
