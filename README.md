# rampr -- Ramp Ripening 


### Overview 
This project provides a framework for simulating a model for [particle
ripening](http://journal.frontiersin.org/article/10.3389/fphy.2014.00018/full)
with non-conserved overall volume. The framework consists of a `C` library
`rampr` and a corresponding Python package `pyrampr`. 

The `C` library `rampr` provides performance-sensitive routines that can be
used to simulate the ripening model. It can solve the dynamics with both an
Euler and an embedded Runge-Kutta-Cash-Karp integrator.
The Python package `pyrampr` employs `ctypes` to access the functionality
of `rampr`. It provides higher level routines that make conducting diverse
simulations comfortable.

### Installation
In order to use `rampr` or `pyrampr`, the shared library `librampr.so` must
be compiled first. On modern unix-like systems, a simple `make` should be
sufficient.

If the compilation was sucessfull, `pyrampr` should be usable via `import
pyrampr` in a python shell. Note that `pyrampr` depends on the packages
`numpy` and `scipy`, which must be installed independently.  When trying to
import `pyrampr` from other directories than the root of this repository,
the root folder should be added to `PYTHONPATH` first.

### Usage
The most important part of `pyrampr` is the class `Simulation`, defined in
`pyrampr/core.py`. Instances of this class must be created by providing an
initial radii distribution, i.e. an array of positive float-values. Each
value corresponds to the radius of one particle. The second ingredient
needed to construct a `Simulation` is a value for the parameter `k`, which
determines the strength of the volume growth.
```python
# Create an initial distribution of radii
# In this case, the 1000 radii are uniformly distributed between 1 and 3
radii = np.linspace(1, 3, 1000)

# Create a simulation object with growth parameter 3.5
sim = Simulation(radii, 3.5)
```
This simulation can now be evolved in time, until either a certain
time, particle number, or volume is reached.
```python
# Make a single evolution step with the Runge-Kutta-Cash-Karp algorithm
sim.evolve_rkck()

# Evolve until the intern time is 2.5
sim.evolve_rkck(until_time=2.5)
```

A number of properties of the simulation can be obtained as member
variables:
```python
sim.radii       # The current radii distribution
sim.volumes     # The current volume distribution
sim.nr_living   # Number of droplets that have not evaporated
sim.nr_dead     # Number of droplets that have evaporated
sim.time        # Simulation time
sim.xi          # Volume growth strength
sim.k           # Volume growth parameter
sim.rmean       # Mean radius
```
For a complete list, look at the `@properties` defined in
`pyrampr/core.py`.


