import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from random import randint

def inOrbit(r = 6378, z = 500):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')

    # nesting functions since they will use the same figure, axes setup
    def plotEarth(r = r):   # can use the parameter given in parent function as default

        # creating meshgrid of phi/theta angles for sphere creation
        pi = np.pi
        phi, theta = np.mgrid[0: pi: 100j, 0: 2*pi: 100j]

        # equations of a sphere in coordinate space
        radX = r*np.sin(phi)*np.cos(theta)
        radY = r*np.sin(phi)*np.sin(theta)
        radZ = r*np.cos(phi)

        # plotting earth model
        ax.plot_surface(radX, radY, radZ, alpha = 0.5)

    # creating empty list for point
    xPts = []
    yPts = []
    
    # creating path/list for point (circle)
    # NOTICE: THIS "ORBIT" DOES NOT USE ORBITAL MECHANICS AND IS JUST AN ARBITRARY CIRCLE AROUND THE EARTH!
    orbitRad = r + z
    theta = np.linspace(0, 2*np.pi, 101)
    xList = np.array([orbitRad*np.cos(thetaPt) for thetaPt in theta])
    yList = np.array([orbitRad*np.sin(thetaPt) for thetaPt in theta])

    # animation function (must be within loop())
    def animate(i):
        # appending new points
        xPts.append(xList[i])
        yPts.append(yList[i])

        ax.clear()  # clearing old plot
        plotEarth() # plotting earth model
        ax.plot(xList, yList, 0, color = "blue", label = "Orbit")    # plotting path
        ax.scatter(xPts[-1],yPts[-1], 0, color = "red", label = "Spacecraft")    # plotting new point

        # setting window parameters
        ax.set_xlim([-r, r])
        ax.set_ylim([-r, r])
        ax.set_zlim([-r, r])
        ax.set_box_aspect((1,1,1))

        # plot labels
        ax.legend()
        ax.set_xlabel('X [km]')
        ax.set_ylabel('Y [km]')
        ax.set_zlabel('Z [km]')
        ax.set_title('Model of Earth centered at origin')

        # plot cleanliness (less clutter with label ticks)
        plt.locator_params(axis='x', nbins=5)
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='z', nbins=5)

    ani = FuncAnimation(fig, animate, frames = 100, interval = 50, repeat = True)

    # saving output animation as GIF
    file = r"/Users/gagandeepthapar/Desktop/Projects/interactiveOrbitVisualizer/testAnimationOrbit.gif"
    writergif = PillowWriter(fps = 60)
    ani.save(file, writer = writergif)

    # displaying plot
    plt.show()
    
def loop():
    # instantiating figure
    fig = plt.figure()
    ax = plt.axes(projection = '3d')

    # creating empty list for point
    xPts = []
    yPts = []
    
    # creating path/list for point (circle)
    theta = np.linspace(0, 2*np.pi, 101)
    xList = np.array([np.cos(thetaPt) for thetaPt in theta])
    yList = np.array([np.sin(thetaPt) for thetaPt in theta])

    # animation function (must be within loop())
    def animate(i):
        # appending new points
        xPts.append(xList[i])
        yPts.append(yList[i])

        ax.clear()  # clearing old plot
        ax.plot(xList, yList, 0, label = "Path")    # plotting path
        ax.scatter(xPts[-1],yPts[-1], 0, color = "red", label = "Point")    # plotting new point
        ax.set_xlim([-1, 1]) # setting limits
        ax.set_ylim([-1, 1])
        ax.legend() # showing legend
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Point moving on path")

        plt.locator_params(axis='x', nbins=5)
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='z', nbins=5)

    ani = FuncAnimation(fig, animate, frames = 100, interval = 50, repeat = True)

    def onclick(event, pauseFlag):
        
        if pauseFlag:
            ani.pause()
            print(f"PAUSED!")
        else:
            ani.resume()
            print(f"RESUMED!")

        print(f"PAUSE FLAG: {pauseFlag}")
        
        pauseFlag = not pauseFlag

    cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, True))
    print(cid)

    # fig.canvas.mpl_connect('button_press_event', toggle_pause(paused))


    # displaying plot
    plt.show()
    
if __name__ == '__main__':
    loop()
    # inOrbit()