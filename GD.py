from Integration import *
from Base import *
from numpy import *
from Flux import *

class GD:
    
    #Galerkin discontinu
    def __init__(self,ordre,omega):
        self.ordre = ordre
        self.omega = omega
        self.N = len(self.omega)
        self.flux = Flux(self.ordre,self.omega)
        self.integ = Integration(self.ordre,self.omega)
        self.base = Base(self.ordre)

        self.A = arange(4,dtype=float)
        self.A = self.A.reshape(2,2)
        self.A[0][0] = 0
        self.A[1][0] = 1
        self.A[0][1] = 1
        self.A[1][1] = 0
        
        self.M = arange(4,dtype=float)
        self.M = self.M.reshape(2,2)
        self.M[0][0] = 1
        self.M[1][0] = 0
        self.M[0][1] = 0
        self.M[1][1] = 1
        
        return

    def MatriceM(self,maille):
        M = arange(self.ordre*self.ordre,dtype=float)
        M = M.reshape(self.ordre,self.ordre)
        
        for i in range(self.ordre):
            for j in range(self.ordre):
                M[i][j] = self.integ.integphilphim(maille,i,j)
        
        return M
    
    def MatriceK(self,maille):
        K = arange(self.ordre*self.ordre,dtype=float)
        K = K.reshape(self.ordre,self.ordre)

        for i in range(self.ordre):
            for j in range(self.ordre):
                K[i][j] = self.integ.integphilphimp(maille,i,j)
        
        return K
    
    def GDm(self,maille,Wmm1,Wm,Wmp1):
        M = self.MatriceM(maille)
        #print("M = ",M)        
        K = self.MatriceK(maille)
        #print("K = ",K)
        
        #F = self.flux.MatriceF(maille,Wmm1,Wm,Wmp1)
    
        #print("Wm",maille,Wm)
        
        B = -dot(K,Wm)
        #print("B",B)
        
        dtW = numpy.linalg.solve(M,B)
        
        return dtW
    
    def GD(self,W):
        i=0
        while i<(self.N-1):  
            w0 = arange(self.ordre*2,dtype=float)
            w0 = w0.reshape(self.ordre,2)
            w1 = arange(self.ordre*2,dtype=float)
            w1 = w1.reshape(self.ordre,2)
            w2 = arange(self.ordre*2,dtype=float)
            w2 = w2.reshape(self.ordre,2)
            
            if i==0:
                for m in range(self.ordre):
                    w0[m][0] = W[(i)*self.ordre+m][0]
                    w0[m][1] = W[(i)*self.ordre+m][1]
                    w1[m][0] = W[i*self.ordre+m][0]
                    w1[m][1] = W[i*self.ordre+m][1]
                    w2[m][0] = W[(i+1)*self.ordre+m][0]
                    w2[m][1] = W[(i+1)*self.ordre+m][1]
            else:
                if i==self.N-2:
                    for m in range(self.ordre):
                        w0[m][0] = W[(i-1)*self.ordre+m][0]
                        w0[m][1] = W[(i-1)*self.ordre+m][1]
                        w1[m][0] = W[i*self.ordre+m][0]
                        w1[m][1] = W[i*self.ordre+m][1]
                        w2[m][0] = W[(i)*self.ordre+m][0]
                        w2[m][1] = W[(i)*self.ordre+m][1]
                else:
                    for m in range(self.ordre):
                        w0[m][0] = W[(i-1)*self.ordre+m][0]
                        w0[m][1] = W[(i-1)*self.ordre+m][1]
                        w1[m][0] = W[i*self.ordre+m][0]
                        w1[m][1] = W[i*self.ordre+m][1]
                        w2[m][0] = W[(i+1)*self.ordre+m][0]
                        w2[m][1] = W[(i+1)*self.ordre+m][1]
                                    
            w = self.GDm(i,w0,w1,w2)
            
            for m in range(self.ordre):
                W[i*self.ordre+m][0] = w[m][0] 
                W[i*self.ordre+m][1] = w[m][1]
            i+=1
        return W