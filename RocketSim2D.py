# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math as m
import launchSimUtilities as lsu


class Rocket(object):
    """Describes the specifications of a rocket. This includes staging, engine specs, payloads and (eventually)
    aerodynamic behaviors."""

    def __init__(self, stages, initial_position):
        try:
            initial_position.shape
        except Exception:
            raise Exception('In object \'Rocket\', \'initial_position\' must be numpy.array of shape (2,)')
        self.stages = stages
        self.position = initial_position  # TODO: investigate array shape manipulations
        self.velocity = np.array([0, lsu.earthGroundSpeed([lsu.earthRadius, 0, 0])])  # TODO: this should sense initial velocity based on initial_position (longitude)
        self.pitch = 90

    def total_mass(self):
        """Returns the total mass of the rocket"""

        mass = 0
        for ii in self.stages:
            if ii.attached:
                mass += ii.dry_mass + ii.fuel_mass
        return mass

    def weight(self):
        """Returns the weight force vector of the whole rocket"""

        xhat = self.position/np.linalg.norm(self.position)
        return -xhat*self.total_mass()*lsu.g(self.position)


class Stage(object):
    """Describes a single stage of a rocket. Takes engine objects as inputs and
    is passed to a rocket object."""

    def __init__(self, name, cross_section, drag_coefficient, engine, num_engines, dry_mass, fuel_mass):
        super(Stage, self).__init__()
        self.name = name
        self.cross_section = cross_section
        self.drag_coefficient = drag_coefficient
        self.engine = engine
        self.num_engines = num_engines
        self.dry_mass = float(dry_mass)
        self.fuel_mass = float(fuel_mass)
        self.active = False
        self.attached = True
        self.final_stage = False

    def thrust(self, R):
        """Returns the thrust provided by a certain engine stage"""

        thrust_magnitude = self.num_engines*lsu.g0*self.engine.Isp(R)*self.engine.massrate  # Thrust magnitude
        # thrustVector = np.array([[0.70710678], [0.70710678]])
        # thrust_vector = np.array([1, 0])
        # return thrust_vector*thrust_magnitude
        return thrust_magnitude


    def drag(self, R, V):  # TODO: remove the portion of drag from theta velocity corresponding to earthly rotation
        """returns the drag force exerted on a certain stage."""

        if lsu.altitude(R) > 200000:
            return np.array([0, 0])
        else:
            if np.linalg.norm(V) == 0:
                return np.array([0, 0])
            else:
                vhat = V/np.linalg.norm(V)  # correct for speed of rotation of Earth
                V = np.linalg.norm(V)
                return (-0.5*lsu.density(R)*(V**2)*self.drag_coefficient*self.cross_section)*vhat


class Engine(object):
    """Describes the specifications of an engine including Isp, thrust, fuel_mass
    flow, fuel_mass type and other things that I will eventually get to."""

    def __init__(self, Isp_vac, Isp_sea, massrate):
        super(Engine, self).__init__()
        self.Isp_vac = float(Isp_vac)
        self.Isp_sea = float(Isp_sea)
        self.massrate = float(massrate)
        # self.thrust = float(lsu.g0*Isp*massrate)

    def Isp(self, R):
        if lsu.altitude(R) < 80000:
            return self.Isp_sea + (1/lsu.pressure(0))*(lsu.pressure(0)-lsu.pressure(R))*(self.Isp_vac-self.Isp_sea)
        else:
            return self.Isp_vac


class Telemetry(object):
    """
    Tracks the flight of a rocket and all the parameters of the flight.
    Telemetry should be a class different than rocket so that I can pass this to the eventual "real" simulator as an ideal trajectory.
    takes in a rocket object. makes blank position, velocity, etc. vectors for that rocket. tracks engine cutoffs, stage separations
    """
    def __init__(self, vehicle, coast_time, dt):
        super(Telemetry, self).__init__()
        self.vehicle = vehicle
        self.total_time = 0
        for activeStage in vehicle.stages:
            if activeStage.fuel_mass != 0:
                dm = activeStage.engine.massrate*activeStage.num_engines*dt  # kg/sec
                self.total_time += (activeStage.fuel_mass*dt)/dm
            else:
                self.total_time += coast_time
        self.position = np.empty((2, self.total_time/dt))
        self.velocity = np.empty((2, self.total_time/dt))
        self.acceleration = np.empty((2, self.total_time/dt))
        self.force = np.empty((2, self.total_time/dt))
        # TODO: add an orientation variable. LVLH, FlightPath angle
        # use velocity vector as a reference in LVLH. horizontal plane, local vertical vector, velocity vector
        # orientation should be given in two angles with respect to... something...

    def telemetry_plot(self):
        t = np.linspace(0, self.total_time, len(self.position[0]))
        plt.figure(1)
        plt.subplot(221)
        plt.title('position')
        plt.plot(self.position[0], self.position[1], 'g')
        # plt.plot(meco[1],meco[0]-lsu.earthRadius,'x')
        plt.subplot(222)
        plt.title('velocity')
        plt.plot(self.velocity[0], self.velocity[1])
        plt.subplot(223)
        # plt.title('acceleration')
        # plt.plot(self.acceleration[0], self.acceleration[1])
        plt.title('force')
        plt.plot(self.force[0], self.force[1])
        plt.subplot(224)
        plt.title('elevation')
        plt.plot(t, np.linalg.norm(self.position, axis=0) - lsu.earthRadius)
        plt.show()

# TODO: create control class to control thrust vector, PITCH KICK
# takes telemetry as an argument
# class trajectory(object):
#     """Takes a rocket object and final orbit parameters as arguments, then
#     computes an ideal trajectory to reach that orbit. trajectory.flight
#     simulates the flight of the rocket guided by the ideal trajectory."""
#     def __init__(self, periapsis, apoapsis, inclination, rocket, flightTime, dt):
#         super(trajectory, self).__init__()
#         self.periapsis = float(periapsis)
#         self.apoapsis = float(apoapsis)
#         self.inclination = m.radians(float(inclination))
#         self.rocket = rocket
#         self.flightTime = flightTime
#         self.dt = dt
class Control(object):
    """
    Control contains functions and parameters for steering the rocket.
    """
    # There will be two functions.
    # One contains a list of instructions to control thrust vectoring duing an ideal flight
    # Another will sense current position relative to an ideal trajectory and make corrections. Used in "real" simulations

    def __init__(self, vehicle, telemetry):
        self.vehicle = vehicle
        self.telemetry = telemetry

    def maneuver_list(self):
        if lsu.altitude(self.telemetry.position) > 30000:
            self.vehicle.pitch = m.radians(0)
        elif lsu.altitude((self.telemetry.position)) <= 30000:
            self.vehicle.pitch = m.radians(90)
        frame_angle = m.atan(self.vehicle.position[0]/self.vehicle.position[1])
        theta = self.vehicle.pitch - frame_angle
        return np.array([m.cos(theta), -m.sin(theta)])

def maneuver_list(vehicle):
    if lsu.altitude(vehicle.position) > 25000:
        vehicle.pitch = m.radians(0)
    elif lsu.altitude(vehicle.position) <= 25000:
        vehicle.pitch = m.radians(90)
    try:
        frame_angle = m.atan(vehicle.position[0]/vehicle.position[1])
    except RuntimeWarning:  # TODO hack on this...
        frame_angle = 0
    theta = vehicle.pitch - frame_angle
    return np.array([m.cos(theta), -m.sin(theta)])


# class  control(object):
#     thrust vector
#     Thrust direction vector
#     thrust magnitude (throttle)
#     orientation
#     self.pitch = m.radians(90) # orientation of vehicle to local horizontal.


def flight(vehicle, telemetry, coast_time, dt):
    """Simulates the flight of a vehicle."""
    a = np.array([[0], [0]])
    F = np.array([[0], [0]])
    total_time = 0
    posn = vehicle.position  # define initial position, velocity, etc.
    velo = vehicle.velocity
    accl = np.array([[0], [0]])
    forc = vehicle.weight()
    for active_stage in vehicle.stages:
        active_stage.active = True
        print active_stage.name, ' start position: ', vehicle.position
        print active_stage.name, ' start velocity: ', vehicle.velocity
        print active_stage.name, ' force at start: ', F
        print active_stage.name, ' mass at start: ', vehicle.total_mass()
        if active_stage.fuel_mass != 0:
            dm = active_stage.engine.massrate*active_stage.num_engines*dt  # kg/sec
            T = (active_stage.fuel_mass*dt)/dm
            run_time_counter = T/dt
        else:
            run_time_counter = coast_time/dt
            dm = 0
        for ii in range(int(total_time), int(total_time+run_time_counter)):
            if active_stage.fuel_mass > 0:
                F = (maneuver_list(vehicle)*active_stage.thrust(vehicle.position)) + active_stage.drag(vehicle.position, vehicle.velocity) + vehicle.weight()
                # TODO: consider having a vectoring term from control object, while thrust variable from Stage is a constant.
                # Vectoring term could also adjust throttle.
                active_stage.fuel_mass -= dm
            elif active_stage.fuel_mass <= 0:
                F = active_stage.drag(vehicle.position, vehicle.velocity) + vehicle.weight()
            else:
                print 'something went wrong'
            a = F/vehicle.total_mass()
            vehicle.velocity = vehicle.velocity + (a*dt)
            vehicle.position = vehicle.position + (vehicle.velocity*dt)  # + (1/2)*a*(dt**2) # the acceleration term. may not be necessary.
            telemetry.position[..., ii] = vehicle.position
            telemetry.velocity[..., ii] = vehicle.velocity
            telemetry.acceleration[..., ii] = a
            telemetry.force[..., ii] = F
        total_time += run_time_counter
        print active_stage.name, ' end position: ', vehicle.position
        print active_stage.name, ' end velocity: ', vehicle.velocity
        print active_stage.name, ' force at end: ', F
        print active_stage.name, ' mass at end: ', vehicle.total_mass()
        if active_stage.final_stage:
            active_stage.active = False
        else:
            active_stage.active = False
            active_stage.attached = False
    return [posn, velo, accl, forc]  # , totalTime]
