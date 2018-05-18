'''
Program wide conventions
R denotes the position vector of the spacecraft
V denotes the velocity vector of the spacecraft
default units will be m, kg, and seconds
'''
import math as m
import numpy as np

## Constants
earthRadius = 6371008.0  # m
earthMass = 5.9721986e24  # kilograms
G = 6.67384e-11  # (N*m^2)/(kg^2) Universal Gravity constant
earthMu = G * earthMass
g0 = earthMu / (earthRadius ** 2)
g = lambda R: earthMu / (np.linalg.norm(R) ** 2)
earthGroundSpeed = lambda R: ((2 * m.pi) / 86400) * m.sin(m.acos(R[2] / np.linalg.norm(R))) * earthRadius

# index) altitude,
# 0) Temperature
# 1) Acceleration of gravity
# 2) Absolute pressure
# 3) Density
# 4) Dynamic viscosity
atmosphere = {-1000: [21.50, 9.810, 11.39, 13.47, 1.821],
              0: [15.00, 9.807, 10.13, 12.25, 1.789],
              1000: [8.50, 9.804, 8.988, 11.12, 1.758],
              2000: [2.00, 9.801, 7.950, 10.07, 1.726],
              3000: [-4.49, 9.797, 7.012, 9.093, 1.694],
              4000: [-10.98, 9.794, 6.166, 8.194, 1.661],
              5000: [-17.47, 9.791, 5.405, 7.364, 1.628],
              6000: [-23.96, 9.788, 4.722, 6.601, 1.595],
              7000: [-30.45, 9.785, 4.111, 5.900, 1.561],
              8000: [-36.94, 9.782, 3.565, 5.258, 1.527],
              9000: [-43.42, 9.779, 3.080, 4.671, 1.493],
              10000: [-49.90, 9.776, 2.650, 4.135, 1.458],
              15000: [-56.50, 9.761, 1.211, 1.948, 1.422],
              20000: [-56.50, 9.745, 0.5529, 0.8891, 1.422],
              25000: [-51.60, 9.730, 0.2549, 0.4008, 1.448],
              30000: [-46.64, 9.715, 0.1197, 0.1841, 1.475],
              40000: [-22.80, 9.684, 0.0287, 0.03996, 1.601],
              50000: [-2.5, 9.654, 0.007978, 0.01027, 1.704],
              60000: [-26.13, 9.624, 0.002196, 0.003097, 1.584],
              70000: [-53.57, 9.594, 0.00052, 0.0008283, 1.438],
              80000: [-74.51, 9.564, 0.00011, 0.0001846, 1.321]}


def rotation(theta, axis):
    """This function creates a rotation transformation matrix about the specified axis, X, Y, or Z through an angle theta"""

    if axis == 'X' or axis == 'x':
        rotationMatrix = np.matrix([[1, 0, 0],
                                    [0, m.cos(theta), -m.sin(theta)],
                                    [0, m.sin(theta), m.cos(theta)]])
    elif axis == 'Y' or axis == 'y':
        rotationMatrix = np.matrix([[m.cos(theta), 0, -m.sin(theta)],
                                    [0, 1, 0],
                                    [m.sin(theta), 0, m.cos(theta)]])
    elif axis == 'Z' or axis == 'z':
        rotationMatrix = np.matrix([[m.cos(theta), m.sin(theta), 0],
                                    [-m.sin(theta), m.cos(theta), 0],
                                    [0, 0, 1]])
    else:
        print('Not a valid rotation axis. Or other error')
    return rotationMatrix


def pressure(R):
    alt = min(atmosphere.keys(),
              key=lambda i: abs(i - altitude(R)))  # get the altitude index in atmostphere closest to actual altitude
    return atmosphere[alt][2]


def density(R):
    """Return the density of the air at current altitude."""

    alt = min(atmosphere.keys(),
              key=lambda i: abs(i - altitude(R)))  # get the altitude index in atmostphere closest to actual altitude
    return atmosphere[alt][3]


def altitude(R):
    return np.linalg.norm(R) - earthRadius


def cart2polar(R):
    x = R[0]
    y = R[1]
    # z = R[2]  #3D
    # r = m.sqrt( x**2 + y**2 + z**2 )  # 3D
    r = m.sqrt( x**2 + y**2 )
    # phi = m.acos(z/r)  #3D
    theta = m.atan(x/y)
    # return np.array([r, theta, phi])  # 3D
    return np.array([r, theta])

def geo2lvlh(vector, vehicle_position):
    x = vehicle_position[0]
    y = vehicle_position[1]
    z = vehicle_position[2]  # 3D
    # phi = m.atan(x/z)  # 3D
    theta = m.atan(y/x)
    # return rotation(phi, 'y') * rotation(theta, 'z') * vector  # 3D
    h = rotation(theta, 'z') * vector.reshape(3, 1)
    return h.reshape(3,) - vehicle_position

def lvlh2geo(vector, vehicle_position):
    x = vehicle_position[0]
    y = vehicle_position[1]
    z = vehicle_position[2]  # 3D
    # phi = -m.atan(x/z)  # 3D
    theta = -m.atan(y/x)
    # return rotation(theta, 'z') * rotation(phi, 'y') * vector  # 3D
    h = rotation(theta, 'z') * vector.reshape(3, 1)
    return h.reshape(3,) + vehicle_position