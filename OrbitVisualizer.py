import numpy as np
from numpy.linalg import solve
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# class creation for arbitrary central body
class CentralBody:
    def __init__(self, name, radius, Mu = None, mass = None):
        self.name = name
        self.rad = radius
        if(Mu != None):
            self.mu = Mu
        if(mass != None):
            G = 6.6743015*(10**-11)
            self.mu = G*mass

    def __repr__(self):
        return (f"{self.name}: Radius {self.rad} km; Mu {self.mu} km^3*s^-2")

# test function - do not use for createOrbit
def twoBodySol(R, V):

    def twoBody(t, state):
        R = np.array([state[0], state[1], state[2]])
        V = np.array([state[3], state[4], state[5]])

        rad = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
        
        ax = -1*398600*R[0] / rad**3
        ay = -1*398600*R[1] / rad**3
        az = -1*398600*R[2] / rad**3

        dState = np.array([V[0], V[1], V[2], ax, ay, az])

        return dState
    
    state = np.array([R[0], R[1], R[2],
                        V[0], V[1], V[2]])
                        
    return solve_ivp(twoBody, (0, 2*3600), state, method = 'RK23', t_eval = np.linspace(0, 2*3600, 1001))

# main driver code to plot orbit and spacecraft orbiting the body
def createOrbit(RInitial, VInitial, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
    # instantiate figure
    fig = plt.figure()
    ax = plt.axes(projection = '3d')

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

    hBar = np.cross(RInitial, VInitial) # specific angular momentum of orbit

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
    ani = FuncAnimation(fig, animate, frames = 1000, interval = np.ceil(1000*60/period), repeat = True)
    
    plt.show()

if __name__ == '__main__':
    R = np.array([6878, 0, 0])
    V = np.array([0, np.sqrt(398600/R[0]), 0])
    createOrbit(R, V)