import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.patheffects as path_effects
from _helperFuncs import *

class interactiveAnimation:
    def __init__(self, RInitial, VInitial, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
        self.R = RInitial
        self.V = VInitial
        self.body = centralBody

        self.computeBody()
        self.computePath()

        self.fig = plt.figure()
        self.ax = plt.axes(projection = '3d')

        # flags
        self.pauseFlag = False
        self.legendFlag = True
        self.bodyFlag = True

        return

    def __repr__(self):
        return f"\nCustom Interactive Animation Class with R-Initial: {self.R} km; V-Initial: {self.V} km/s\n"

    def computePath(self):
        def twoBody(t, state):
            curR = np.array([state[0], state[1], state[2]])
            curV = np.array([state[3], state[4], state[5]])

            rad = np.linalg.norm(curR)

            ax = -1*self.body.mu*curR[0]/(rad**3)
            ay = -1*self.body.mu*curR[1]/(rad**3)
            az = -1*self.body.mu*curR[2]/(rad**3)

            dState = np.array([curV[0], curV[1], curV[2], ax, ay, az])

            return dState
        
        state = np.array([self.R[0], self.R[1], self.R[2], self.V[0], self.V[1], self.V[2]])

        self.semiMajor = -self.body.mu/(np.linalg.norm(self.V)**2 - (2*self.body.mu/np.linalg.norm(self.R))) # equation for semi-major axis
        self.period = 2*np.pi*(self.semiMajor**1.5)/(np.sqrt(self.body.mu)) # equation for period of orbit

        sol = solve_ivp(twoBody, (0, np.ceil(self.period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(self.period), 1000))
        
        self.orbitState = np.array([sol.t,
                            sol.y[0], sol.y[1], sol.y[2],
                            sol.y[3], sol.y[4], sol.y[5]])

        return self.orbitState

    def computeBody(self):
        pi = np.pi
        phi, theta = np.mgrid[0:pi:100j, 0:2*pi:100j]

        self.bodyX = self.body.rad*np.sin(phi)*np.cos(theta)
        self.bodyY = self.body.rad*np.sin(phi)*np.sin(theta)
        self.bodyZ = self.body.rad*np.cos(phi)

    def plotSystem(self):
        def animate(i):
            self.ax.clear()
            if self.bodyFlag:
                self.ax.plot_surface(self.bodyX, self.bodyY, self.bodyZ, alpha = 0.5)
            self.ax.plot(self.orbitState[1], self.orbitState[2], self.orbitState[3], color = 'blue', label = 'Orbit')
            self.ax.scatter(self.orbitState[1][i], self.orbitState[2][i], self.orbitState[3][i], color = 'red', label = 'Spacecraft')

            # set window boundaries
            self.ax.set_xlim([-1*self.semiMajor, self.semiMajor])
            self.ax.set_ylim([-1*self.semiMajor, self.semiMajor])
            self.ax.set_zlim([-1*self.semiMajor, self.semiMajor])

            # plot labels
            if self.legendFlag:
                self.ax.legend()    
            self.ax.set_xlabel('X [km]')
            self.ax.set_ylabel('Y [km]')
            self.ax.set_zlabel('Z [km]')
            self.ax.set_title('Spacecraft Orbiting Earth')
            
            # plot cleanliness (less clutter with label ticks)
            self.ax.set_box_aspect((1, 1, 1))
            plt.locator_params(axis='x', nbins=5)
            plt.locator_params(axis='y', nbins=5)
            plt.locator_params(axis='z', nbins=5)

        self.ani = FuncAnimation(self.fig, animate, frames = 1000, interval = np.ceil(1000*60/self.period), repeat = True)
        self.fig.canvas.mpl_connect('key_press_event', self.userInputs)
        
        plt.show()

    def userInputs(self, event):
        # pause animation
        if event.key == 'p':
            if self.pauseFlag:
                self.ani.resume()
            else:
                self.ani.pause()
            
            self.pauseFlag = not self.pauseFlag
        
        # show legend
        if event.key == 'l':
            self.legendFlag = not self.legendFlag
        
        if event.key == 'b':
            self.bodyFlag = not self.bodyFlag
        


if __name__ == '__main__':
    R = np.array([6878, 0, 0])
    V = np.array([0, np.sqrt(398600/(2*R[0])), np.sqrt(398600/(2*R[0]))])
    anim = interactiveAnimation(R, V)
    anim.plotSystem()
    # print(anim.orbitState[0])

    # A = np.array([1, 2, 3, 4])
    # B = np.array([5, 6,7, 8])
    # C = np.append([A, B, A, B, A, B])
    # print(C)
    