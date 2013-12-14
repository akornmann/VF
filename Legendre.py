import numpy

class Legendre:
    def __init__(self,ordre):
        self.ordre = ordre
        self.w1 = numpy.array([2.0],dtype=float)
        self.w2 = numpy.array([1.0,1.0],dtype=float)
        self.w3 = numpy.array([5.0/9.0,8.0/9.0,5.0/9.0],dtype=float)
        self.wa = numpy.array([self.w1,self.w2,self.w3])
        
        self.x1 = numpy.array([0.0],dtype=float)
        self.x2 = numpy.array([-1.0/numpy.sqrt(3.0),1.0/numpy.sqrt(3.0)],dtype=float)
        self.x3 = numpy.array([-numpy.sqrt(3.0/5.0),0.0,numpy.sqrt(3.0/5.0)],dtype=float)
        self.xa = numpy.array([self.x1,self.x2,self.x3])
        
    def w(self):
        return self.wa[self.ordre-1]
    
    def x(self):
        return self.xa[self.ordre-1]