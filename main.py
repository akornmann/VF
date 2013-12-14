from GD import *
from Integration import *
from Legendre import *
from RK2 import *
from animation import *

ordre=3
dx=0.01
dt=dx/(2*(ordre+1))
tmax=1.0
a = -2.0
b = 10.0


omega = numpy.arange(a,b,dx)
N = len(omega)

g = Integration(ordre,omega)
l = Legendre(ordre)
xa = l.x()
W = arange(ordre*2*(N-1),dtype=float)
W = W.reshape((N-1)*ordre,2)

#Initialisation
x=[]
l=0
while l<N-1:
    for i in range(ordre):
        xv=g.G(xa[i],l+1)
        x.append(xv)

        W[l*ordre+i][0] = exp(-xv*xv)
        W[l*ordre+i][1] = exp(-xv*xv)
    l+=1

r = RK2(ordre,W,omega)

t=0
W = r.Wp1(W, dt)


anim = Animation(x,-1.5,1.5)
while t<tmax:
    #W = r.Wp1(W, dt)
    t+=dt
    l=[]
    
    i=0

    while i<ordre*(N-1):
        l.append(W[i].flat[1])
        i+=1
    anim.plot(l)

        