from numpy import polynomial

class Legendre:
    def __init__(self,ordre):
        (self.x, self.w) = polynomial.legendre.leggauss(ordre)
 
    def weight(self):
        return self.w
    
    def coord(self):
        return self.x