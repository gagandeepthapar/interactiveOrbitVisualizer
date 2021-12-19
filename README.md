# **Interactive Orbit Visualizer**
Plot an orbit around the Earth. Control your spacecraft. See how your orbit changes!

## **Motivation**
Using Python libraries to create a model of an arbitrary orbit and using user input to change that orbit via impulse burns.

I wanted to showcase some orbital mechanics and software development skills in the same project without the looming goal of a class grade.

I strutured the repository with different files to showcase a style of development I typically walk through - by creating small, sometimes irrelevant, use cases of a unique function (i.e., 3D plotting) to understand how a feature works before incorporating it in a larger setting where the requirements and implementation are more concrete.

## **File Structure**
### **testPlot.py:**</br>
``testPlot.py`` explored matplotlib's ability to plot in 3-Dimensional space. The script models the Earth as a sphere (with radius **6378 km**) and a basic orbit (modeled as a circle with radius **6878km**).</br></br>

The ``testPlot.png`` file shows the output of the test code.</br></br>

### **testAnimation.py:**
``testAnimation.py`` explored matplotlib's ability to animate a 3-Dimensional plot. The script created two instances: a simple loop and a simple orbit. The loop script modeled a circular path with a single point traversing it over the length of the animation. The orbit script modeled the Earth and a basic orbit with a spacecraft propagating over that orbit. It should be noted that the orbit was *not* constructed using Kepler's Equations of Motion, but as a simple circle with a radius greater than that of the Earth.</br></br>

The ``testAnimationLoop.gif`` is the output animation of the loop script while ``testAnimationOrbit.gif`` is the output animation of the orbit script.</br></br>

### **OrbitVisualizer.py**</br>
``OrbitVisualizer.py`` is the culmination of the discoveries of the previous test files. This script takes either the **state vector** (position and velocity) of a given spacecraft or the **Classical Orbital Elements** to create an acccurate representation of the spacecraft's orbit around the Earth. The orbit is created using a Runge Kutta (2,3) solver to integrate across Kepler's Equations of Motion. </br></br>
Users can also interact with the model using the arrow keys to perform impulse maneuvers. Each maneuver represents a 500 m/s (0.5 km/s) impulse burn in a specified direction and see how their orbit is altered.</br></br>

The ``OrbitVisualizer.gif`` file is the creation of a simple 500 km-altitude circular orbit with 0 deg inclination. It should also be noted that the propagation of the spacecraft is 1 sec : 1 min. That is, 1 second of the animation represents 1 minute of the spacecraft's actual propagation.

## **Notable Libraries**
* NumPy
* matplotlib
* SciPy