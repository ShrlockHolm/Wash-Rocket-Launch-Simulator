import numpy as np
from math import *
import launchSimUtilities as lsu
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

steps = 100
t = np.linspace(0, 2*pi, steps) # time vector

class orbit:
	'''The parameters defining an ideal keplerian orbit. For now we ignore the
	position of the orbiting body.'''
	def __init__(self, periapsis, apoapsis, raan, inc, argPeri):
		#super(orbit, self).__init__()
		self.periapsis = float(periapsis)		# km
		self.apoapsis = float(apoapsis)			# km
		self.raan = float(raan)					# Right Ascention of Ascending Node, degrees
		self.inc = float(inc)					# orbit Iinclination, degrees
		self.argPeri = float(argPeri)			# Argument of Periapsis, degrees

		self.a = (self.periapsis+self.apoapsis)/2;		# Semimajor axis, km
		self.e = (self.a-self.periapsis)/self.a;		# eccentricity
		self.b = self.a * sqrt(1-self.e**2);			# Semiminor axis,km

	def plot(self):
		x = self.a*np.cos(t) - (self.a*self.e);
		y = self.b*np.sin(t);
		z = 0*t;
		R = np.array([ [x],
					   [y],
					   [z] ]) # This gives the orbit as if it were in 2D

		# Now we need to rotate that same orbit to a 3D position using THIS:
		R = lsu.rotation(radians(-self.raan),'z') * lsu.rotation(radians(-self.inc),'x')* lsu.rotation(radians(-self.argPeri),'z') * R
		R = R.A # make the matrix R an array so that it plots properly

		X = R[0,...]
		Y = R[1,...]
		Z = R[2,...]

		#Here we actually plot the orbit
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.plot3D(X,Y,Z)
		ax.plot3D([0],[0],[0],'x')

		plt.show()
