from numpy import *
from Base import *

class Flux:
    def __init__(self,ordre,omega):
        self.ordre = ordre
        self.omega = omega
        self.N = len(omega)
        self.Ap = matrix('0.5 0.5;0.5 0.5')
        self.Am = matrix('-0.5 0.5;0.5 -0.5')
        self.b = Base(self.ordre)
        
    def MatriceF(self,maille,Wkm1,Wk,Wkp1):
        A = arange(self.ordre*self.ordre,dtype=float)
        A = A.reshape(self.ordre,self.ordre)
        B = arange(self.ordre*self.ordre,dtype=float)
        B = A.reshape(self.ordre,self.ordre)

        if(maille==0):
            for i in range(self.ordre):
                for j in range(self.ordre):
                    A[i][j] = self.b.phi(self.omega[maille],i)*self.b.phi(self.omega[maille],j)
                    B[i][j] = 0 #a revoir    
                      
            return dot(B,Wkm1)+dot(A,Wk)
        else:
            if(maille==self.N-1):
                for i in range(self.ordre):
                    for j in range(self.ordre):
                        A[i][j] = self.b.phi(self.omega[maille],i)*self.b.phi(self.omega[maille],j)
                        B[i][j] = self.b.phi(self.omega[maille-1],i)*self.b.phi(self.omega[maille],j) #a revoir        
                return dot(B,Wkm1)+dot(A,Wk)
            else:
                for i in range(self.ordre):
                    for j in range(self.ordre):
                        A[i][j] = self.b.phi(self.omega[maille],i)*self.b.phi(self.omega[maille],j)
                        B[i][j] = self.b.phi(self.omega[maille-1],i)*self.b.phi(self.omega[maille],j) #a revoir        
                return dot(B,Wkm1)+dot(A,Wk)