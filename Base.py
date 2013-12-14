from numpy import *
from Legendre import *

class Base:
    
    #ordre : degre des polynomes de la base
    def __init__(self,ordre):
        self.ordre = ordre
        
        leg = Legendre(self.ordre)
        self.w = leg.w()
        self.x = leg.x()
        
        return
    
    #Produit des xi-xj (i!=j)
    def pden(self, d):
        res = 1
        
        for k in range(self.ordre):
            if(k!=d):
                res*=self.x[d]-self.x[k]
        
        return res
    
    #base des phi
    def phi(self,x,d):
        res = 1
        
        for k in range(self.ordre):
            if(k!=d):
                res*=x-self.x[k]
        
        res/=self.pden(d)
        
        return res
    
    #derive de phi
    def phiprime(self,x,d):
        k=0
        l=0
        res=0
        tmp=1
        for k in range(self.ordre):
            if k!=d:
                for l in range(self.ordre):
                    if l!=k:
                        tmp*=x-self.x[l]
                    l+=1
                res+=tmp
                tmp=1
                l=0
            k+=1

        k=0
        tmp=1
        for j in range(self.ordre):
            if k!=j:
                tmp*= (self.x[j]-self.x[k])
            j+=1
        
        res/=tmp
        
        return res