# Object tracking with Bayesian Filtering systems: Kalman filter and Particle filter

As part of a course in Statistical Signal Processing.

## Description
A Kalman filter and Particle Filter implementation for tracking an object in 1) a Gaussian linear model (GaussianModel_ObjectTracking.m) and in the 2) FitzHugh-Nagumo model for the action potential in a cell membrane (ActionPotential_ObjectTracking.m).

![fitz](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/fitz_model.png)

(Example of the FitzHugh-Nagumo model)



### Linear model

For a Gaussian linear model, the Kalman filter provides the optimal solution. Despite the Particle filter not providing the optimal solution, due to its approximations, the solution is still accurate enough to track the object in the plane. There are still some peaks in the error curves and some simulations where the estimate deviates significantly from the real trajectory (due to the stochastic nature of the algorithm).

![trajectory](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/Linear_trajectory.png)

![velocities](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/Linear_velocities.png)



### FitzHugh-Nagumo model

For a non-linear model, only the Particle Filter provides a solution.

![traject-fitz](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/Fitz_trajectory.png)


![traject-fitz2d](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/Fitz_2d.png)


![traject-fitz_estimate](https://github.com/CesarCaramazana/BayesianFiltering/blob/main/Figures/Fitz_estimate.png)

## Conclusions
In this laboratory we have analyzed the role of the parameters for a linear-Gaussian system when
estimating the state with both a Kalman filter and a particle filter. We discovered that, for a simplistic
model, such as an object moving in a 2D space, the difference between the Kalman solution (optimal)
and the particle filter is not that significant.
For the FitzHugh-Nagumo model, only the particle filter could be implemented. We found out the
effect that the number of particles have on the error –inversly proportional–, and that the form of
the observations is critical for the proper functioning of the filter, as with non-linear observations the
problem could not be solved at all.
