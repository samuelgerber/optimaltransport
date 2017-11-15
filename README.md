#Optimal Transport

![scale-4.png](https://bitbucket.org/repo/XyGX46/images/333242785-scale-4.png)
![scale-6.png](https://bitbucket.org/repo/XyGX46/images/661701334-scale-6.png)
![scale-11.png](https://bitbucket.org/repo/XyGX46/images/104944428-scale-11.png)

This code is designed for fast multiscale optimal transport but also permits single scale optimal transport computations and more:

* Optimal transport solved with Lemon, CPLEX, GPLK or MOSEK
* Fast approximate transport using a mutliscale strategy 
* Sinkhorn optimal transport approach
* Multiscale version of sinkhorn transport
* Rpackage for sinkhorn transport (both multiscale and single scale, depends on gmra R package)
* Rpackage for optimal transport (both multiscale and single scale, requires any of CPLEX, GPLK or MOSEK, depends on gmra r package )

## Requirements

The build process is using CMake.

R packages:

* gmra https://github.com/suppechasper/gmra

C/C++ Libraries:

* Eigen
* For some of the comandline  version CPLEX, MOSEK or GLPK (currently setup for CPLEX but few modifications are required to change it to MOSEK or GLPK, edit the .cxx driver files in commandline accordingly to use the different solver header files and instantiate the solver you want to use. Then edit CMake files to point to the library specific header files and libraries. )

## R Package mop Installation

Note: The package just switched to the lemon network simplex optimizer by
default. It has not been thoroughly tested. The most extensively tested setup
is with the CPLEX network optimizer.

The R package requires the gmra package:

* gmra https://github.com/suppechasper/gmra

The package can be installed with

* make R_mop_install (after running cmake)

Prepackaged versions are available from

* https://bitbucket.org/suppechasper/optimaltransport/downloads/

### Using other linear programming libraries

The default setup uses lemon for solving linear progragraming. If you would
like to use another linear programming the pakage supports CPLEX, MOSEK and
GLPK. 

To setup for a a specific library edit src/Makevars and src/mop_config.h in the
Rpackage/mop directory before buidling the package.

To install the package requires that all the libraries you are using are on R's linker path.
