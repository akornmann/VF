from pylab import ion,plot,axis,draw


class Animation:
    def __init__(self, x, ymin, ymax):
        self.x = x
        self.ymin = ymin
        self.ymax = ymax
        ion()
        self.line = []
        l, = plot(x,x)
        self.line.append(l)
        l2, = plot(x,x)
        self.line.append(l2)
        axis([x[0], x[-1], ymin, ymax])
   
    def plot(self, y):
        for i in range(len(self.line)):
            self.line[i].set_xdata(self.x)
            self.line[i].set_ydata(y[i])

        draw()