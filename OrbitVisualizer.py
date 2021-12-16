import numpy as np
from matplotlib import pyplot as plt

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

def plotCentralBody(centralBody):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')

    pi = np.pi
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]

    # equations of a sphere in coordinate geometry
    radX = centralBody.rad*np.sin(phi)*np.cos(theta)
    radY = centralBody.rad*np.sin(phi)*np.sin(theta)
    radZ = centralBody.rad*np.cos(phi)

    # plotting
    ax.plot_surface(radX, radY, radZ, alpha = 0.5)
    ax.set_xlim([-centralBody.rad, centralBody.rad])
    ax.set_ylim([-centralBody.rad, centralBody.rad])
    ax.set_zlim([-centralBody.rad, centralBody.rad])

    ax.set_xlabel("X [km]")
    ax.set_ylabel("Y [km]")
    ax.set_zlabel("Z [km]")

    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='z', nbins=5)

    ax.set_box_aspect((1,1,1))
    plt.show()

def createOrbit(R, V, centralBody = CentralBody("Earth", 6378, 3.986*(10**5))):
    # create orbit path from R, V vectors

    # plot orbit path and central body
    plt = plotCentralBody(centralBody)
    plt.plot([-centralBody.rad, centralBody.rad], [-centralBody.rad, centralBody.rad], [-centralBody.rad, centralBody.rad])

    plt.show()

if __name__ == '__main__':
    createOrbit(1,2)