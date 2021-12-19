import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from _helperFuncs import *

class interAnimation:
    # initialization
    def __init__(self, RInitial, VInitial, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
        # setting pause flag
        self.paused = False

        # instantiate figure
        self.fig = plt.figure()
        self.ax = plt.axes(projection = '3d')

        fig = self.fig
        ax = self.ax

        # FUNCTIONS:
        # create orbit path from R, V vectors
        def twoBody(t, state):
            # unpacking state vector
            R = np.array([state[0], state[1], state[2]])
            V = np.array([state[3], state[4], state[5]])

            # definition of radius 
            rad = np.linalg.norm(R)

            # Orbital Kinematic Equations
            ax = -1*centralBody.mu*R[0] / (rad**3)
            ay = -1*centralBody.mu*R[1] / (rad**3)
            az = -1*centralBody.mu*R[2] / (rad**3)

            dState = np.array([V[0], V[1], V[2], ax, ay, az])

            return dState
        
        # instantiating state vector
        state = np.array([R[0], R[1], R[2],
                            V[0], V[1], V[2]])

        a = -centralBody.mu/(np.linalg.norm(VInitial)**2 - (2*centralBody.mu/np.linalg.norm(RInitial))) # equation for semi-major axis
        period = 2*np.pi*(a**1.5)/(np.sqrt(centralBody.mu)) # equation for period of orbit

        # running ODESolver to get radius, velocity vectors
        # solving for tSpan of [0, period] with 1000 points in between
        sol = solve_ivp(twoBody, (0, np.ceil(period)), state,method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1000))

        # position vectors across the period
        spaceX = sol.y[0]
        spaceY = sol.y[1]
        spaceZ = sol.y[2]

        # plot orbit path and central body
        def plotCentralBody(centralBody = centralBody):
            # creating meshgrid of phi/theta angles for sphere creation
            pi = np.pi
            phi, theta = np.mgrid[0: pi: 100j, 0: 2*pi: 100j]

            # equations of a sphere in coordinate space
            radX = centralBody.rad*np.sin(phi)*np.cos(theta)
            radY = centralBody.rad*np.sin(phi)*np.sin(theta)
            radZ = centralBody.rad*np.cos(phi)

            # plotting earth model
            ax.plot_surface(radX, radY, radZ, alpha = 0.5)

            # set aspect ratio
            ax.set_box_aspect((1,1,1))

        # PLOTTING
        def animate(i):
            ax.clear()  # clear old plot
            plotCentralBody()   # plot earth
            ax.plot(spaceX, spaceY, spaceZ, color = 'blue', label = 'Orbit')    # plot orbit path
            ax.scatter(spaceX[i], spaceY[i], spaceZ[i], color = 'red', label = 'Spacecraft')    # plot point on orbit
        
            # setting window parameters
            r = np.linalg.norm(RInitial)
            ax.set_xlim([-r, r])
            ax.set_ylim([-r, r])
            ax.set_zlim([-r, r])
            ax.set_box_aspect((1,1,1))

            # plot labels
            ax.legend()
            ax.set_xlabel('X [km]')
            ax.set_ylabel('Y [km]')
            ax.set_zlabel('Z [km]')
            ax.set_title('Spacecraft Orbiting Earth\n1 sec = 1 min')

            # plot cleanliness (less clutter with label ticks)
            plt.locator_params(axis='x', nbins=5)
            plt.locator_params(axis='y', nbins=5)
            plt.locator_params(axis='z', nbins=5)

        # animating plot
        # NOTICE: The resulting animation has a 1:1 sec:min ratio; that is, 1 min in reality is 1 sec in the animation
        self.ani = FuncAnimation(fig, animate, frames = 1000, interval = np.ceil(1000*60/period), repeat = True)

        self.fig.canvas.mpl_connect('key_press_event', self.togglePause)
    
    # custom toString
    def __repr__(self):
        return f"Custom Interactive Animation Class"
        
    # custom pause function
    def togglePause(self, event):
        # only pauses when 'p' is pressed
        if event.key == 'p':
            if self.paused:
                self.ani.resume()
            else:
                self.ani.pause()
            
            self.paused = not self.paused

    

if __name__ == '__main__':
    R = np.array([6878, 0, 0])
    V = np.array([0, np.sqrt((398600/R[0])/2), np.sqrt((398600/R[0])/2)]) # inc = 45deg; z = 500km

    iAni = interAnimation(R, V)
    plt.show()