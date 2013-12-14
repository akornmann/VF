from numpy  import *
from GD import *

class RK2:
    
    def __init__(self,ordre,W0,omega):
        self.ordre=ordre
        self.omega = omega
        self.W0 = W0
        self.GD = GD(self.ordre,self.omega)
    
    def dt(self,W):
        return self.GD.GD(W)
    
    def Wp1(self,W, deltat):
        dtW = self.dt(W)
        We = self.sum(W,self.scal(dtW,0.5*deltat))
        dtWe = self.dt(We)
        Wp1 = self.sum(W,self.scal(dtWe,deltat))

        return Wp1
    
    def sum(self,A,B):
        if(len(A)!=len(B)) :
            print("Sizing error")
        C=[]
        for i in range(len(A)):
            C.append(A[i]+B[i])
        return C
    
    def scal(self,A,l):
        C=[]
        for i in range(len(A)):
            C.append(l*A[i])
        return C