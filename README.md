Incremental fluids Java

This is a port to Java, resorting to Apache Commons Math, with some minor changes of the C++ code from https://github.com/tunabrain/incremental-fluids

The port is not a strait forward, as instead of single files I separated some code in common among them into classes:
-The solids are now individual classes common to the fluid simulators classes.
-The linear interpolators were gathered in a single class stored into a maths package, integrators could possibly be also stored there something to consider.
-I made an exhaustive effort to maintain the files additions among solvers by keeping even comments among Fluid solvers files methods and variables.

I do not make promises but there may appear a 2D versions a few more solvers (LBM, SWE, Stable fluids, Practical Fluids, and SPH). 
Also, there may appear a 3D versions of some of these codes someday.
Finally, parallel versions either in CUDA/OpenCL/GLSL. 

I would gladly accept suggestions/contributions to make the code more clearer or to add interesting features.
