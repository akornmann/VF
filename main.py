from numpy import zeros,exp,sin,pi,sqrt
from GD import GD
from RK2 import RK2
from Animation import Animation


'''
def identity(W):
    return W

dt=0.001
rk2 = RK2(dt,identity)
rk2.unit()
'''
'''
from Base import Base
b = Base(3)
b.unit()
'''

#param galerkin
N = 20
ordre = 4
a=0.
b=2.

#param rk2
dt=((b-a)/N)/(2*(ordre+1))
Tmax=20.

#param animation
ymin = -1.0
ymax = 1.0

#init galerkin
galerkin = GD(ordre,a,b,N)
#Array 3D contenant la solution
W = zeros((2,N,ordre))

#init animation
x = galerkin.getX()
anim = Animation(x,ymin,ymax)

#init rk2
rk2 = RK2(dt,galerkin.dW)

#avance temporelle
while(rk2.t<Tmax):
    
    def lim0():
        we = zeros((2,1))
        #we[:,0] = [10*exp(-100*(rk2.t-0.2)*(rk2.t-0.2)),0.]
        we[:,0] = [sin(5*rk2.t/pi),0.]
        return we
    
    def limN():
        we = zeros((2,1))
        we[:,0] = [0.,0.]
        return we
    
    galerkin.Wexact0 = lim0
    galerkin.WexactN = limN
    
    W = rk2.Wp1(W)
    y = galerkin.getY(W)
    
    anim.plot(y)

from Integration import Integration
integ = Integration(ordre)
norme = 0
for i in range(N):
    for j in range(ordre):
        norme+=(W[0,i,j]-sin(5*x[i*(ordre)+j]/pi))*(W[0,i,j]-sin(5*x[i*(ordre)+j]/pi))
        
norme = sqrt(norme)

print norme
    