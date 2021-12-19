import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

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

# class creation for orbit around arbitrary body
class Orbit:
    def __init__(self, stateVector = None, coesList = None, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
        if stateVector.any() == None and coesList.any() == None:
            print(f"\nPlease input either a State Vector [R, V] or a list of Classical Orbital Elements with the following format:\n")
            print("Angular Momentum [kg3*m-2]\nEccentricity\nInclination [deg]\nRAAN [deg]\nArgument of Perigee [deg]\nTrue Anomaly [deg]\n")
        else:
            if(coesList == None):
                self.R = stateVector[0]
                self.V = stateVector[1]

                self.h, self.ecc, self.inc, self.raan, self.arg, self.ta, self.a = RV2coes(self.R, self.V, centralBody)

            else:
                self.h, self.ecc, self.inc, self.raan, self.arg, self.ta = coesList

                self.R, self.V = coes2RV(coesList)
        
        self.body = centralBody
    
    def __repr__(self):
        ra = self.h ** 2 / self.body.mu * (1/(1 + self.ecc))
        rp = self.h ** 2 / self.body.mu * (1/(1 - self.ecc))
        
        a = ra + rp

        beta = np.arccos(-1*self.ecc)
        b = a * ((1-self.ecc**2)/(1 + self.ecc * np.cos(beta)))

        return f"{round(a, 3)} km x {round(b,3)} km orbit with inclination of {self.inc} deg"

    def createPath(self):
        sol = twoBodySol(self.R, self.V, self.body)
        rX = sol.y[0]
        rY = sol.y[1]
        rZ = sol.y[2]
        vX = sol.y[3]
        vY = sol.y[4]
        vZ = sol.y[5]

        self.rPath = np.array([rX, rY, rZ])
        self.vPath = np.array([vX, vY, vZ])

        return np.array([self.rPath, self.vPath])

# ODE solver for Kepler's Equations of Motion
def twoBodySol(R, V, centralBody = CentralBody("Earth", 6378, Mu = 398600)):

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

    a = -centralBody.mu/(np.linalg.norm([V[0], V[1], V[2]])**2 - (2*centralBody.mu/np.linalg.norm([R[0], R[1], R[2]]))) # equation for semi-major axis
    period = 2*np.pi*(a**1.5)/(np.sqrt(centralBody.mu)) # equation for period of orbit

                        
    return solve_ivp(twoBody, (0, np.ceil(period)), state, method = 'RK23', t_eval = np.linspace(0, np.ceil(period), 1001))

# returning state vector from classical orbital elements
def coes2RV(angularMomentum, eccentricity, inclination, raan, argumentOfPerigee, trueAnomaly, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
    # Explanation of parameters
        # inc - [deg] inclination of orbit
        # raan - [deg] right ascenscion of ascending node
        # ecc - [~] eccentricity of orbit
        # arg - [deg] argument of perigee
        # h - [km3s2] specific angular momentum of orbit
        # ta - [deg] true anomaly of spacecraft on orbit
        # a - [km] semi major axis of orbit
        # M - [deg] mean anomaly (different from true anomaly)
    # renaming parameters for use
    h = angularMomentum
    ecc = eccentricity
    raan = np.deg2rad(raan)
    inc = np.deg2rad(inclination)
    arg = np.deg2rad(argumentOfPerigee)
    ta = np.deg2rad(trueAnomaly)

    periRConst = (h**2/centralBody.mu) * (1/(1 + ecc * np.cos(ta)))
    periVConst = centralBody.mu/h

    periRBar = periRConst * np.array([np.cos(ta), np.sin(ta), 0])
    periVBar = periVConst * np.array([-1*np.sin(ta), ecc + np.cos(ta), 0])

    R = 0
    V = 0

    Qa = np.array([-np.sin(raan)*np.cos(inc)*np.sin(arg) + np.cos(raan)*np.cos(arg), 
                    -np.sin(raan)*np.cos(inc)*np.cos(arg) - np.cos(raan)*np.sin(arg),
                        np.sin(raan)*np.sin(inc)])

    Qb = np.array([np.cos(raan)*np.cos(inc)*np.sin(arg) + np.sin(raan)*np.cos(arg),
                    np.cos(raan)*np.cos(inc)*np.cos(arg) - np.sin(raan)*np.sin(arg),
                    -np.cos(raan)*np.sin(inc)])
    
    Qc = np.array([np.sin(inc)*np.sin(arg), np.sin(inc)*np.cos(arg), np.cos(inc)])

    Q = np.array([Qa, Qb, Qc]).reshape(3,3)

    R = np.matmul(Q, periRBar)
    V = np.matmul(Q, periVBar)

    return R, V

# returning classical orbital elements from state vector
def RV2coes(R, V, centralBody = CentralBody("Earth", 6378, Mu = 398600)):
    # Explanation of parameters
        # R - [km, km, km] position vector of spacecraft at some time, t
        # V - [km/s, km/s, km/s] velocity vector of spacecraft at time t
    # Explanation of variables
        # h - angularMomentum [km3s-2]
        # inc - inclination [deg]
        # ecc - eccentricity
        # raan - right ascenscion of ascending node [deg]
        # arg - argumentOfPerigee [deg]
        # ta - trueAnomaly [deg]
        # a - semiMajorAxis [km]]

    rad = np.linalg.norm(R)
    vel = np.linalg.norm(V)
    
    radVel = np.dot(R, V)/rad
    
    hBar = np.cross(R, V)
    h = np.linalg.norm(hBar)

    inc = np.rad2deg(np.arccos(hBar[2]/h))

    nodeBar = np.cross(np.array([0, 0, 1]), hBar)
    node = np.linalg.norm(nodeBar)

    if node == 0:
        raan = 0
    else:
        raan = np.rad2deg(np.arccos(nodeBar[0]/node))
        if nodeBar[1] < 0:
            raan = 360-raan
    
    eccBar = 1/centralBody.mu * ((vel**2 - centralBody.mu/rad) * R - rad*radVel*V)
    ecc = np.linalg.norm(eccBar)

    if node == 0:
        arg = 0
    else:
        arg = np.rad2deg(np.arccos(np.dot(nodeBar/node , eccBar/ecc)))
        if eccBar[2] < 0:
            arg = 360-arg
    
    ta = np.rad2deg(np.arccos(np.dot(eccBar/ecc , R/rad)))
    if radVel < 0:
        ta = 360-ta
    
    rP = (h**2 / centralBody.mu) * (1/(1 + ecc))
    rA = (h**2 / centralBody.mu) * (1/(1 - ecc))
    
    a = 0.5*(rP + rA)

    # return [angularMomentum [km3s-2], eccentricity, inclination [deg], raan [deg], argumentOfPerigee [deg], trueAnomaly [deg], semiMajorAxis [km]]
    return np.array([h, ecc, inc, raan, arg, ta, a])

if __name__ == '__main__':
    R = np.array([6878, 0, 0])
    V = np.array([0, np.sqrt(398600/6878), 0])
    o = Orbit(stateVector=np.array([R,V]))