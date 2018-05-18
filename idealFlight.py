'''
Program wide conventions
R denotes the position vector of the spacecraft
V denotes the velocity vector of the spacecraft
default units will be m, kg, and seconds
'''
import numpy as np
import math as m
import launchSimUtilities as lsu
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class rocket(object):
    """Describes the specifications of a rocket. This includes staging, engine
    specs, payloads and (eventually) aerodynamic behaviors."""
    def __init__(self, stages, latitude, longitude):
		self.stages = stages
		self.V = np.array([surface_vel(latitude),0])
		self.lat = latitude
		self.vsuf = np.array([surface_vel(latitude),0])
		self.pitch = math.pi/2

class stage(object):
    """Describes a single stage of a rocket. Takes engine objects as inputs and
    is passed to a rocket object."""
	def __init__(self, name, crossSection, dragCoefficient, engine, numEngines, drymass, fuel):
        super(stage, self).__init__()
		self.name = name
		self.crossSection = crossSection
		self.dragCoefficient = dragCoefficient
		self.engine = engine
		self.numEngines = numEngines
		self.drymass = drymass
		self.fuel = fuel
		self.active = False
		self.attached = True
    def thrust(R):
        return self.numEngines*lsu.g0*self.engine.Isp(R)*self.engine.massrate

class engine(object):
    """Describes the specifications of an engine including Isp, thrust, fuel
    flow, fuel type and other things that I will eventually get to."""
    def __init__(self, Isp_vac, Isp_sea,massrate):
        super(engine, self).__init__()
        self.Isp_vac = float(Isp_vac)
        self.Isp_sea = float(Isp_sea)
        self.massrate = float(massrate)
        self.thrust = float(lsu.g0*Isp*massrate)
    def Isp(R):
    	if altitude(R) < 80000:
    		return self.Isp_sea + (1/lsu.pressure(0))*(lsu.pressure(0)-lsu.pressure(R))*(self.Isp_vac-self.Isp_sea)
    	else:
    		return Isp_vac


class trajectory(object):
    """Takes a rocket object and final orbit parameters as arguments, then
    computes an ideal trajectory to reach that orbit. trajectory.flight
    simulates the flight of the rocket guided by the ideal trajectory."""
    def __init__(self, periapsis, apoapsis, inclination, rocket, flightTime, dt):
        super(trajectory, self).__init__()
        self.periapsis = float(periapsis)
        self.apoapsis = float(apoapsis)
        self.inclination = m.radians(float(inclination))
        self.rocket = rocket
        self.flightTime = flightTime
        self.dt = dt

    def flight(arg):


        for i, stage in enumerate(self.stages):
			if stage.active and stage.fuel < 0:
				print "Stage cutoff!"
				stage.active = False
				stage.attached = False
				try: self.stages[i+1].active = True
				except: pass
