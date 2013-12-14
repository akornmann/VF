#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This class can be used to display results of 1D solver.
"""


import pylab as pl

class Animation:
    """2D animated plot.
    """
    
    # public interface
    
    def __init__(self, x, ymin, ymax, *args):
        """ x : array of abscissae
            ymin, ymax : min and max of ordinates
        """
      
        self.x = x
        
        self.ymin = ymin
        self.ymax = ymax


        pl.ion()
        
        self.line, = pl.plot(x[0], x[-1], *args)
        
        pl.axis([x[0], x[-1], ymin, ymax])
    
    
    
    def plot(self, y):
      """Update line with given y values.
      """
      self.line.set_xdata(self.x)
      self.line.set_ydata(y)
      pl.draw()
