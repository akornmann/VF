from Base import *

class Integration:
    
    #ordre : degre des polynomes de la base
    #maille : numero de la maille courante
    #omega
    def __init__(self,ordre,omega):
        self.ordre = ordre
        self.omega = omega
        
        leg = Legendre(self.ordre)
        self.w = leg.w()
        self.x = leg.x()
        
        self.b = Base(self.ordre)
        return
    
    def integphilphim(self,maille,l,m):
        res = 0
        
        for k in range(self.ordre):
            x = self.x[k]
            fx = self.b.phi(x,l)*self.b.phi(x,m)
            res+=self.G(fx,maille)*self.w[k]

        res = res*2*(self.omega[3]-self.omega[2]) #norme de la maille
        
        return res
    
    def integphilphimp(self,maille,l,m):
        res = 0
        
        for k in range(self.ordre):
            x = self.x[k]
            fx = self.b.phi(x,l)*self.b.phiprime(x,m)
            res+=self.G(fx,maille)*self.w[k]

        res = res*2*(self.omega[3]-self.omega[2]) #norme de la maille
        
        return res
    
    def G(self,x,maille):
        return (self.omega[maille-1]+self.omega[maille])/2+x*(self.omega[maille]-self.omega[maille-1])/2
    
    def Gm1(self,x,maille):
        return (x-(self.omega[maille-1]+self.omega[maille])/2)*2/(self.omega[maille]-self.omega[maille-1])