# Object tracking with Bayesian Filtering systems: Kalman filter and Particle filter

As part of a course in Statistical Signal Processing.

## Description
A Kalman filter and Particle Filter implementation for tracking an object in 1) a Gaussian linear model (GaussianModel_ObjectTracking.m) and in the 2) FitzHugh-Nagumo model for the action potential in a cell membrane (ActionPotential_ObjectTracking.m).



## Conclusions
In this laboratory we have analyzed the role of the parameters for a linear-Gaussian system when
estimating the state with both a Kalman filter and a particle filter. We discovered that, for a simplistic
model, such as an object moving in a 2D space, the difference between the Kalman solution (optimal)
and the particle filter is not that significant.
For the FitzHugh-Nagumo model, only the particle filter could be implemented. We found out the
effect that the number of particles have on the error –inversly proportional–, and that the form of
the observations is critical for the proper functioning of the filter, as with non-linear observations the
problem could not be solved at all.
