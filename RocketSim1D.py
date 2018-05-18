from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math as m
import launchSimUtilities as lsu

class rocket(object):
	"""Describes the specifications of a rocket. This includes staging, engine
	specs, payloads and (eventually) aerodynamic behaviors."""
	def __init__(self, stages, initialPosition):
		# self.activeStage = stages[0]
		# self.nextStages = stages[1:]
		self.stages = stages
		# self.V = np.array([lsu.earthGroundSpeed(latitude),0,0])
		self.velocity = 0
		self.position = initialPosition
		self.pitch = m.pi/2
	def totalMass(self):
		'''Returns the total mass of the rocket'''
		mass = 0
		for ii in self.stages:
			if ii.attached == True:
				mass += ii.dryMass + ii.fuelMass
		return mass
	def weight(self):
		'''Returns the weight force vector of the whole rocket'''
		xhat = self.position/np.linalg.norm(self.position)
		return -xhat*self.totalMass()*lsu.g(self.position)

class stage(object):
	"""Describes a single stage of a rocket. Takes engine objects as inputs and
	is passed to a rocket object."""
	def __init__(self, name, crossSection, dragCoefficient, engine, numEngines, dryMass, fuelMass):
		super(stage, self).__init__()
		self.name = name
		self.crossSection = crossSection
		self.dragCoefficient = dragCoefficient
		self.engine = engine
		self.numEngines = numEngines
		self.dryMass = float(dryMass)
		self.fuelMass = float(fuelMass)
		self.active = False
		self.attached = True
	def thrust(self, R):
		'''Returns the thrust provided by a certain engine stage'''
		return self.numEngines*lsu.g0*self.engine.Isp(R)*self.engine.massrate
	def drag(self, R, V):
		'''returns the drag force exerted on a certain stage.'''
		if lsu.altitude(R) > 200000:
			return 0
		else:
			if np.linalg.norm(V) == 0:
				return 0
			else:
				vhat = V/np.linalg.norm(V) # reduce the velocity vector to its components to accurately point the drag vector
				V = np.linalg.norm(V)
				return (-0.5*lsu.density(R)*(V**2)*self.dragCoefficient*self.crossSection)*vhat

class engine(object):
	"""Describes the specifications of an engine including Isp, thrust, fuelMass
	flow, fuelMass type and other things that I will eventually get to."""
	def __init__(self, Isp_vac, Isp_sea,massrate):
		super(engine, self).__init__()
		self.Isp_vac = float(Isp_vac)
		self.Isp_sea = float(Isp_sea)
		self.massrate = float(massrate)
		#self.thrust = float(lsu.g0*Isp*massrate)
	def Isp(self, R):
		if lsu.altitude(R) < 80000:
			return self.Isp_sea + (1/lsu.pressure(0))*(lsu.pressure(0)-lsu.pressure(R))*(self.Isp_vac-self.Isp_sea)
		else:
			return self.Isp_vac

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

def flight(Vehicle, coastTime, dt):
	'''Simulates the flight of a vehicle. '''
	a = 0
	F = 0
	totalTime = 0
	posn = np.array([])
	velo = np.array([])
	accl = np.array([])
	forc = np.array([])
	meco = [0,0]
	for activeStage in Vehicle.stages:
		activeStage.active = True
		# print activeStage.name,' start position: ', Vehicle.position-lsu.earthRadius
		# print activeStage.name,' start velocity: ', Vehicle.velocity
		# print activeStage.name,' force at start: ', F
		# print activeStage.name,' mass at start: ', Vehicle.totalMass()
		if activeStage.fuelMass != 0:
			dm = activeStage.engine.massrate*activeStage.numEngines*dt # kg/sec
			T = (activeStage.fuelMass*dt)/dm
			totalTime += T
			runTime = np.linspace(0,T,T/dt)
		else:
			runTime = np.linspace(0,coastTime,coastTime/dt)
			totalTime += coastTime
		for t in runTime:
			if activeStage.fuelMass>0:
				F = activeStage.thrust(Vehicle.position) + activeStage.drag(Vehicle.position,Vehicle.velocity) + Vehicle.weight()
				activeStage.fuelMass -= dm
			elif activeStage.fuelMass<=0:
				F = activeStage.drag(Vehicle.position,Vehicle.velocity) + Vehicle.weight()
			else:
				print 'something went wrong'
			a = F/(activeStage.fuelMass+activeStage.dryMass)
			Vehicle.velocity += a*dt
			Vehicle.position += Vehicle.velocity*dt + (1/2)*a*(dt**2) # the acceleration term. may not be necessary.
			posn = np.append(posn,Vehicle.position)
			velo = np.append(velo,Vehicle.velocity)
			accl = np.append(accl,a)
			forc = np.append(forc,F)
		# print activeStage.name,' end position: ', Vehicle.position-lsu.earthRadius
		# print activeStage.name,' end velocity: ', Vehicle.velocity
		# print activeStage.name,' force at end: ', F
		# print activeStage.name,' mass at end: ', Vehicle.totalMass()
		activeStage.active = False
		activeStage.attached = False
	return [posn, velo, accl, forc, totalTime]

def telemetryPlot(telemetry):
	posn = telemetry[0]
	velo = telemetry[1]
	accl = telemetry[2]
	forc = telemetry[3]
	totalTime = telemetry[4]
	T = np.linspace(0,totalTime,len(posn))
	plt.figure(1)
	plt.subplot(221)
	plt.plot(T,posn-lsu.earthRadius,'g')
	# plt.plot(meco[1],meco[0]-lsu.earthRadius,'x')
	plt.subplot(222)
	plt.plot(T,velo,'b')
	plt.subplot(223)
	plt.plot(T,accl,'r')
	plt.subplot(224)
	plt.plot(T,forc,'y')
	plt.show()
