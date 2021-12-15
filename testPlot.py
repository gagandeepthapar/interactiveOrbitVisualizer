import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plotCentralBody(r = 6378, zoom = 2000):
    # initiate figure
    fig = plt.figure()
    ax = plt.axes(projection = '3d')

    # create grid of phi, theta values to define a sphere in space
    pi = np.pi
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]

    # equations of a sphere in coordinate geometry
    radX = r*np.sin(phi)*np.cos(theta)
    radY = r*np.sin(phi)*np.sin(theta)
    radZ = r*np.cos(phi)

    # creating trivial orbit
    alt = 1000
    xline = np.linspace(-r - alt, r + alt, (int)((r+alt)/10))
    ylinePOS = np.array([np.sqrt((r+alt)**2 - eleX**2) for eleX in xline])
    ylineNEG = np.array([-1*eleY for eleY in ylinePOS])
    zline = np.array([0 for eleX in xline])

    # plotting
    ax.plot_surface(radX, radY, radZ, alpha = 0.5)
    ax.scatter(0, 0, r, color = "red", label = "North Pole")
    ax.scatter(0, 0, -r, color = "blue", label = "South Pole")
    ax.plot(xline, ylinePOS, zline, color = 'orange')
    ax.plot(xline, ylineNEG, zline, color = 'orange', label = 'Trivial Orbit')

    # setting window limits
    ax.set_xlim([-r - zoom, r + zoom])
    ax.set_ylim([-r - zoom, r + zoom])
    ax.set_zlim([-r - zoom, r + zoom])
    ax.set_box_aspect((1,1,1))

    # defining window properties
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')
    ax.legend()
    ax.set_title('Model of Earth centered at origin')

    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='z', nbins=5)

    # show window
    plt.show()

if __name__ == "__main__":
    plotCentralBody()