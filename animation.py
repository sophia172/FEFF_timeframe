import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import os
import glob
import ffmpeg

current_path = os.getcwd()
os.chdir('/Users/sophia17_2/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/CdS37')

filelist = glob.glob('/Users/sophia17_2/ownCloud/PhD/Matlab figure/CdS/EXAFS/*CdK.dat')

i = 0


fig = plt.figure()
ax = plt.axes(xlim=(0, 20), ylim=(-1, 1))

xdata, ydata = np.array([]), np.array([])
ln, = plt.plot(xdata, ydata, 'k-')

def init():
    ax.set_xlim(0, 20)
    ax.set_ylim(-1, 1)

    return ln,
def update(i):
    print(i)
    data = np.genfromtxt('/Users/sophia17_2/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/CdS37/EXAFSsimulation_CdS{:s}.dat'.format(str(int(i))))[:, (2, 5)]
    xdata = data[:,0]
    ydata = data[:,1] *data[:,0]**2

    ln.set_data(xdata, ydata)
    return ln,

ani = animation.FuncAnimation(fig, update, frames=np.linspace(1,14,14), interval=200)
plt.show()
ani.save('demo.htm', fps=30, extra_args=['-vcodec', 'libx264'])

