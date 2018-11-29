# pulse.py
# Oliver Thomson Brown, Rupert Nash
# 2018-10-31
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import count

import d3q15

# define box size
nx = 10
ny = 10
nz = 50

# set relaxation time
tau = 0.5
amplitude = 0.01

class AniPlotter(object):
    def __init__(self):
        self.lattice = self.init_lattice()
        self.x = np.mgrid[1:nz+1]

        self.fig = plt.figure()
        self.ymax = amplitude
        self.ax = plt.axes(xlim=(0, nz+2), ylim=(-0.5*amplitude, 0.5*amplitude))
        self.line, = self.ax.plot([], [], lw=2)

    @property
    def y(self):
        return self.lattice.rho[nx/2, ny/2, 1:nz+1] - 1.0

    def init_lattice(self):
        lattice = d3q15.Lattice(nx, ny, nz, tau, tau)
        # set boundary conditions and forcing
        lattice.initBoundaryC('periodic')
        lattice.setForceC('none')

        # define initial density
        lattice.rho[:, :, nz/2] += amplitude
        lattice.rho /= lattice.rho[nx/2, ny/2, 1:nz+1].mean()
        lattice.initFromHydroVars()
        return lattice

    def init_plot(self):
        self.line.set_data([], [])
        return self.line,

    def advance_to(self, t):
        self.lattice.step(t - self.lattice.time_step)
        self.lattice.updateHydroVars()

    def get_plot(self, t):
        self.advance_to(t)
        max_amp = 0.5 * self.y.ptp()
        # if max_amp < 0.25 * self.ymax:
        #     self.ymax /= 2.0
        #     self.ax.set_ylim(-self.ymax, self.ymax)
        #     self.ax.figure.canvas.draw()

        self.line.set_data(self.x, self.y)
        return self.line

    def report(self):
        print("Time step ", self.lattice.time_step)
        print("Density:")
        print(self.lattice.rho[nx/2, ny/2, 1:nz+1] - 1.0)
        #print("Velocity:")
        #print(self.lattice.u[nx/2, ny/2, 1:nz+1])

ap = AniPlotter()
#ani = animation.FuncAnimation(ap.fig, ap.get_plot, frames=count(0,5),
#                                  init_func=ap.init_plot, repeat=False, blit=True)
#plt.show()

ap.report()
ap.advance_to(100)
ap.report()
ap.advance_to(200)
ap.report()
