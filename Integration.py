from Legendre import Legendre
from numpy import abs

class Integration:
    def __init__(self,ordre):
        self.ordre = ordre
        
        leg = Legendre(self.ordre)
        self.weigth = leg.weight()
        
        return

    def integral(self,f,a,b):
        L=abs(b-a)
        res = 0.
        for i in range(self.ordre):
            res+=f[i]*self.weigth[i]
            
        res *= (L/2.)
        return res
    
    def Gk(self,x,a,b):
        return ((b-a)*x+(a+b))/2.