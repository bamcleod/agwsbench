import numpy as np
import pylab as pl

nrows = 10000
data = np.zeros((nrows,8))
i=0

for ang0 in [ 90, 180, 270, 0]:
    r = np.random.uniform(6,10,size=nrows) / 60.
    ang = np.random.uniform(ang0-45, ang0+45,size=nrows)
    x = r * np.cos(np.radians(ang))
    y = r * np.sin(np.radians(ang))
    data[:,2*i] = x
    data[:,2*i+1] = y
    pl.plot(x,y,'.')
    i=i+1



np.savetxt('fields.txt', data, fmt='%.4f')
pl.show()
