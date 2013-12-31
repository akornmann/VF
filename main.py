from numpy import zeros,exp
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

N = 50
ordre = 3
a=0.
b=2.

dt=((b-a)/N)/(2*(ordre+1))
Tmax=2.

#Array 3D contenant la solution
W = zeros((2,N,ordre))

galerkin = GD(ordre,a,b,N)

rk2 = RK2(dt,galerkin.dW)

x = galerkin.getX()

anim = Animation(x,-10,10)

while(rk2.t<Tmax):
    def lim0():
        we = zeros((2,1))
        we[:,0] = [15*exp(-100.*(rk2.t-0.3)*(rk2.t-0.3)),0.]
        return we
    galerkin.Wexact0 = lim0  
    
    W = rk2.Wp1(W)
    y = galerkin.getY(W)
    
    anim.plot(y)