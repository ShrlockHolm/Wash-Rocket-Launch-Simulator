import numpy as np
# import math as m
# import matplotlib.pyplot as plt
# import RocketSim2D as rs

class State(object):
    """contains the state vector for a rocket"""

    def __init__(self, position, velocity, time):
        self.position = position
        self.velocity = velocity
        self.time = time

        def x(self):
            return self.position[0]