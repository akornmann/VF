from numpy import ones,zeros,linspace
from Legendre import Legendre

class Base:
    #ordre : degre des polynomes de la base
    def __init__(self,ordre):
        self.ordre = ordre
        
        leg = Legendre(self.ordre)
        self.w = leg.weight()
        self.x = leg.coord()
        return
    
    #base des phi
    def phi(self,x):
        res = ones((self.ordre,1))
        for i in range(self.ordre) :
            for j in range(self.ordre):
                if j!=i:
                    res[i]*=(x-self.x[j])/(self.x[i]-self.x[j]) 
        return res
    
    #derive de phi
    def phip(self, x):
        res = zeros(self.ordre)
        for i in range(self.ordre):
            div=1.
            
            for j in range(self.ordre):
                if j!=i:
                    div*=(self.x[i]-self.x[j])
                    tmp=1.
                    for k in range(self.ordre):
                        if(k!=j and k!=i):
                            tmp*=(x-self.x[k])
                    res[i]+=tmp
            res[i]/=div
        return res

    #test unitaire
    def unit(self):
        import matplotlib.pyplot as plt
        N=500
        
        x=linspace(-1,1,N)
        y=zeros((N,self.ordre))
        for k in range(N):
            y[k] = self.phi(x[k]).T

        plt.plot(x,y)

        plt.title('phi')
        plt.axis([-1,1,-5,5])
        plt.grid()        
        
        plt.show()
        
        
        for k in range(N):
            y[k] = self.phip(x[k]).T

        plt.plot(x,y)

        plt.axis([-1,1,-5,5])
        plt.title('phip') 
        plt.grid()   
        
        plt.show()