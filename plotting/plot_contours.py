from utils import *
from nc_reader import nc_reader
import matplotlib.pyplot as plt

fname = ''
step = 0

ncr = nc_reader()

ncr.open(fname)

ix = 32
iys = 70
iye = 135
izs = 192
ize = 257

v = ncr.get_dataset(step=step, name='y_velocity')
w = ncr.get_dataset(step=step, name='z_velocity')

xg, yg, zg = ncr.get_meshgrid()

ncr.close()


plt.quiver(X=yg[ix, iys:iye, izs:ize],
           Y=zg[ix, iys:iye, izs:ize],
           U=v[ix, iys:iye, izs:ize],
           V=W[ix, iys:iye, izs:ize])

plt.show()


#In line 637 of the paper, I discuss the velocity and vorticity field crossing a front.  I wonder if we could make a zoom showing the behaviour described, maybe as a contour plot with arrows showing the sense of the flow field, and contours showing the vertical vorticity field?  I imagine this will not be that easy, so only attempt it if you think you'd be keen to try!  I think the front seen on the left at t = 63 in fig10 would be interesting to examine; I would take ix = 32 then consider 70 <= iy <= 134 and 192 <= iz <= 256.  This should cross the middle of the front in a yz plane.
