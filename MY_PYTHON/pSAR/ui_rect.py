#!/usr/bin/env python
#
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

class draw_rect(object):
    def __init__(self):
        #
        self.ax      = plt.gca()
        self.x0      = None
        self.y0      = None
        self.x1      = None
        self.y1      = None
        #
        self.polyid  = []
        self.polyflag= 0
        self.pressid = 0
        self.x_index = [2,3]
        self.y_index = [1,2]
        self.poly_x  = np.array([self.x0,self.x0,self.x1,self.x1,self.x0])
        self.poly_y  = np.array([self.y0,self.y1,self.y1,self.y0,self.y0])
        self.polyid_active = 0
        self.npoly         = 0
        self.poly_m        = np.zeros((1,4))
        #self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)

    def on_press(self, event):
        #
        # print(event.button)
        if event.button == 1:
           self.pressid = 1
           self.x0      = event.xdata
           self.y0      = event.ydata
           self.x1      = self.x0
           self.y1      = self.y0
           self.poly_x[:] = self.x0
           self.poly_y[:] = self.y0
           if self.polyflag == 1:
              self.polyid.set_xdata(self.poly_x)
              self.polyid.set_ydata(self.poly_y)
              plt.draw() 
           else:
              self.polyflag  = 1
              self.polyid,   = self.ax.plot(self.poly_x,self.poly_y,'o-r',color='green',lw=3.5, alpha=0.85)
              plt.draw()
        elif event.button == 3:
           self.polyflag = 0
           self.npoly += 1
           if self.npoly == 1:
               self.poly_m = np.array([[np.min(self.poly_x),np.max(self.poly_x),np.min(self.poly_y),np.max(self.poly_y)]])
           else:
               self.poly_m = np.append(self.poly_m,[[np.min(self.poly_x),np.max(self.poly_x),np.min(self.poly_y),np.max(self.poly_y)]],axis=0)
           #
        else:
           # print(dir(plt.Figure))
           # self.ax.show(block=False)
           print(dir(self.ax.figure))
           plt.close()
        #
        
    def on_release(self, event):
        #
        self.pressid = 0 
        return self
        #
    def motion_notify_callback(self,event):
        #
        if self.pressid == 1:
           x, y = event.xdata, event.ydata
           self.poly_x[self.x_index[0]] = x
           self.poly_x[self.x_index[1]] = x
           self.poly_y[self.y_index[0]] = y
           self.poly_y[self.y_index[1]] = y
           self.polyid.set_xdata(self.poly_x)
           self.polyid.set_ydata(self.poly_y)
           plt.draw()       


