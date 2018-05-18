import RocketSim2D as rs
import launchSimUtilities as lsu
import math as m
import numpy as np

Merlin1D = rs.Engine(320, 282, 236)
stage1 = rs.Stage("falcon 9 first stage", m.pi*(3.66/2)**2, 0.3, Merlin1D, 9, 18000, 385000)
stage1.final_stage = True
stage2 = rs.Stage("falcon 9 secnd stage", m.pi*(3.66/2)**2, 0.3, Merlin1D, 1, 4900, 90000)
payload = rs.Stage("payload", m.pi*(3.66/2)**2, 0.3, None, 0, 13500, 0)
falcon9 = rs.Rocket([stage1, stage2, payload], np.array([lsu.earthRadius, 0]))
falcon3 = rs.Rocket([stage1], np.array([lsu.earthRadius, 0]))

falcon9_track = rs.Telemetry(falcon9, 500, 0.05)

rs.flight(falcon9, falcon9_track, 500, 0.05)

falcon9_track.telemetry_plot()
# rs.flight(falcon3, 500, 0.05)
