#### Basic dependencies
#       Metric
#       Forwards & Backward raytracing
#
#
#### Emissivity & absorption Profiles
#    1. Description of spectral density
#    2. Electron Profile
#    3. Velocity Profile
#
#### Radiative Transfer Ray Tracers
#           
#    *Slow light tends to be easier to solve with analytic profiles*
#
#    This is normaly done in 2 steps (For both *Analytic* and Numeric):
#
#    1. Backwards (From observer to BH)
#         Calculates and stores the geodesic embedding
#    2. Forwards (Fom BH to observer)
#        Radiative Transfer usually done with interpolation on stored geodesic embedding
#
#### Test Poblems
#
#    1. 3D non-symmetric flow (Fast Light ϕ integrals)
#          Guassian blob off to the side
#    2. 3D non-symmetric flow (Slow Light ϕ integrals)
#          Gaussian blob on Keplerian orbit
#
#### Ultimate goal (GRMHD simulations)
#    3. 3D (G)RMHD model model ray-trace(mempmap)
#
#
#### Design constraints
#   Keep as much in Julia as possible
#   Threading and GPU friendly (Probably a good idea to not thread along an individual ray)
#   Low allocation on individual raytrace
#
