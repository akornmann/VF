class RK2:
    def __init__(self,dt,F):
        self.dt = dt
        self.F = F
        self.t = 0.
        
    def Wp1(self,W):
        dtW = self.F(W)
        We = W+self.dt/2.*dtW
        dtWe = self.F(We)
        Wp1 = W+self.dt*dtWe
        self.t+=self.dt
        return Wp1
    
    def unit(self):
        from numpy import exp
            
        t=0.    
        W=1.
        for i in range(100):
            W=self.Wp1(W)
            t += self.dt
            Wexact = exp(i*self.dt)
            print("Erreur : ")
            print(W-Wexact)