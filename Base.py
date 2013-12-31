from numpy import ones,zeros,linspace
from Legendre import Legendre

class Base:
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
        from matplotlib.pyplot import plot,title,axis,grid,show
        N=500
        
        x=linspace(-1,1,N)
        y=zeros((N,self.ordre))
        for k in range(N):
            y[k] = self.phi(x[k]).T

        plot(x,y)

        title('phi')
        axis([-1,1,-6,6])
        grid()        
        
        show()
        
        for k in range(N):
            y[k] = self.phip(x[k]).T

        plot(x,y)

        axis([-1,1,-6,6])
        title('phip') 
        grid()   
        
        show()