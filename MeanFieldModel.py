# Preambule
import numpy as np
import matplotlib.pyplot as plt

# Parameters
r = 3.
N = 5.
sigma =1.
f0=0.5
z0=0.5

# Equations
def G(z):
    return (1-r/N)*np.log(z)+(r/2.-1)*np.log(1-z)

def L(f):
    return sigma*np.log(f)+(r-1-sigma)*np.log(1-f)

def H(z,f):
    return G(z)+L(f)
    
def F(z):
    return 1+(r - 1)*z**(N-1)-r/N*(1-z**N)/(1-z)

def df(z,f):
    return -F(z)/(z*(1-z)*(1-z**(N- 1)))

def dz(z,f):
    return (sigma-f*(r-1))/(f*(1-f))

# Time integration
tmax=2
dt=0.0001
fvec=[f0]
zvec=[z0]
tvec=np.linspace(0,tmax,np.int(tmax/dt)+1)
for i in range(0,np.int(tmax/dt)):
    z=zvec[i]
    f=fvec[i]
    
    #z
    k1 = dt*dz(z,f)
    k2 = dt*dz(z+k1/2.,f)
    k3 = dt*dz(z+k2/2.,f)
    k4 = dt*dz(z+k3,f)
    z=z+1./6.*(k1+2.*k2+2.*k3+k4)
    
    #f
    k1 = dt*df(z,f)
    k2 = dt*df(z,f+k1/2.)
    k3 = dt*df(z,f+k2/2.)
    k4 = dt*df(z,f+k3)
    f=f+1./6.*(k1+2.*k2+2.*k3+k4)
    
    
    fvec.append(f)
    zvec.append(z)
    
# Plots
plt.plot(tvec,zvec)
plt.xlabel('Time')
plt.ylabel('Fraction loners of total (z)')
plt.show()
plt.plot(tvec,fvec)
plt.xlabel('Time')
plt.ylabel('Fraction cooperators of participants (f)')
plt.show()