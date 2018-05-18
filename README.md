# WASH v 0.0: Rocket Tracking and Guidance Simulation

This is meant to be a program that will model a rocket launch and perform tracking, as well as guidance on that simulated rocket.

## Program requirements/ steps within the program
1. Build accurate(ish) rocket models
Certain specifications of a real-world rocket will be fed into the program. This includes engine specs, staging, payload, and eventually it will include some aerodynamics.

2. Generate an ideal trajectory
Before the program simulates tracking and guidance it will have to generate an ideal trajectory to reach orbit. I intend to find that ideal trajectory using Monte Carlo methods. The ideal trajectory will be passed to the simulated guidance system as a path to follow.

3. Simulate tracking
Tracking will be performed using simulated tracking systems and Kalman filtering. Real-world tracking systems often have error associated with them. This program will have to simulate the error in addition to performing the same tracking methods.

4. Simulate control
Rocket dynamics will be simulated and guided as directed by the tracking subroutines.

## Project stages
1. Ideal trajectory simulation
First I'll build the ideal trajectory simulator. Most of the rocket dynamics and physics should be worked out in this stage. To find the ideal trajectory I'll have to build a full rocket flight simulator.

2. Tracking
Next, I'll implement the Kalman filter and simulate the tracking sensors by applying Gaussian noise to the ideal trajectory path. In this step the flight will be restricted to the ideal trajectory found in stage 1.

3. Control.
Now I allow the rocket and flight simulator to deviate from the ideal trajectory. I'll apply control vectors to the Kalman filter. The control vectors will include noise to simulate the uncertainty of controls.
