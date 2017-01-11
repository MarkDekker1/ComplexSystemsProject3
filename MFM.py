# Preambule
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Parameters
r = 3.
N = 5.
sigma =1.
x0=0.1
y0=0.6
z0=0.1
f0=x0/(x0+y0)

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
    
# Create Equations
Diction={}
Diction['z']='(sigma-f*(r-1))/(f*(1-f))'
Diction['f']='-1+(r - 1)*z**(N-1)-r/N*(1-z**N)/(1-z)/(z*(1-z)*(1-z**(N- 1)))'
Initials=[z0,f0]

def Equations(state,t=0):
    z,f = state
    list2=[]
    list2.append(eval(Diction['z']))
    list2.append(eval(Diction['f']))
    return list2
    
print 'Computing Integration...'
start = clock()
Results, infodict = integrate.odeint(Equations, Initials, tvec, full_output=True)
print 'done in %.3f seconds!' % (clock()-start)
a,b= Results.T
    
# Plots
plt.plot(tvec,zvec)
plt.xlabel('Time')
plt.ylabel('Fraction loners of total (z)')
plt.show()
plt.plot(tvec,fvec)
plt.xlabel('Time')
plt.ylabel('Fraction cooperators of participants (f)')
plt.show()