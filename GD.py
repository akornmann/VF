from numpy import zeros,linspace,array,dot
from Legendre import Legendre
from Integration import Integration
from Base import Base

class GD:
    #Galerkin discontinu
    def __init__(self,ordre,a,b,N):
        #param
        self.ordre = ordre
        self.a = a
        self.b = b
        self.N = N
        
        #discretisation du segment
        self.ak = linspace(self.a,self.b,self.N+1)
        
        #points et poids d'interpolations Gauss Legendre
        leg = Legendre(self.ordre)
        self.coord = leg.coord()
        self.weight = leg.weight()
        
        #Integration et bijection sur les mailles
        self.integ = Integration(self.ordre)
        #fcts de base
        self.base = Base(self.ordre)
        
        #calcul de phi aux interfaces
        self.phid = self.base.phi(1)
        self.phig = self.base.phi(-1)
        
        #calcul de phip sur les points d'interpolation
        self.dphi = zeros((self.ordre,self.ordre))
        for j in range(self.ordre):
            self.dphi[:,j] = self.base.phip(self.coord[j])
     
    #condition au bord a
    def Wexact0(self):
        we = zeros((2,1))
        return we
    
    #condition au bord b
    def WexactN(self):
        we = zeros((2,1))
        return we
    
    #Calcul de la derivee sptatiale av ec le schema GD
    def dW(self,W):
        dW = zeros((2,self.N,self.ordre))
        for k in range(self.N):
            Lk = self.ak[k+1] - self.ak[k]

            flux = dot(self.flux(k,W),self.phig.T)-dot(self.flux(k+1,W),self.phid.T)            
            for i in range(self.ordre):
                #dot(A,W) est fait directement car integral ne fait que du calcul 1D
                dW[0,k,i]+=self.integ.integral(-W[1,k]*self.dphi[i]*2./Lk, self.ak[k], self.ak[k+1])+flux[0,i]
                dW[1,k,i]+=self.integ.integral(-W[0,k]*self.dphi[i]*2./Lk, self.ak[k], self.ak[k+1])+flux[1,i]
                
                dW[0,k,i] *= 2./Lk/self.weight[i]
                dW[1,k,i] *= 2./Lk/self.weight[i]  
        return dW
    
    #flux a l'interface k                                      
    def flux(self,k,W):
        if(k==0):
            Wk = self.Wexact0()
            Wkp1 = dot(W[:,k],self.phig)
        else:
            if(k==self.N):
                Wk = dot(W[:,k-1],self.phid)
                Wkp1 = self.WexactN()
            else:
                Wk = dot(W[:,k-1],self.phid)
                Wkp1 = dot(W[:,k],self.phig)
        
        return self.F(Wk, Wkp1)
    
    #flux numerique
    def F(self,Wk,Wkp1):
        Amoins = array([[-0.5,-0.5],[-0.5,-0.5]])
        Aplus = array([[0.5,-0.5],[-0.5,0.5]])
        
        return dot(Aplus,Wk) + dot(Amoins,Wkp1)
    
    #formate les valeurs en x pour l'affichage     
    def getX(self):
        formatX=zeros((self.N*self.ordre))
        for i in range(self.N):
            for j in range(self.ordre):
                formatX[i*(self.ordre)+j]= self.integ.Gk(self.coord[j],self.ak[i],self.ak[i+1])

        return formatX
    
    #formate les valeurs en y pour l'affichage
    def getY(self,W):
        formatY=zeros((2,self.N*self.ordre))
        for v in range(2):
            for i in range(self.N):
                for j in range(self.ordre):
                    formatY[v,i*(self.ordre)+j]=W[v,i,j]

        return formatY