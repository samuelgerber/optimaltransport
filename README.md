#Optimal Transport

![scale-4.png](https://bitbucket.org/repo/XyGX46/images/333242785-scale-4.png)
![scale-6.png](https://bitbucket.org/repo/XyGX46/images/661701334-scale-6.png)
![scale-11.png](https://bitbucket.org/repo/XyGX46/images/104944428-scale-11.png)

This code is designed for fast multiscale optimal transport but also permits single scale optimal transport computations and more:

* Optimal transport solved with CPLEX, GPLK or MOSEK
* Fast approximate transport using a mutliscale strategy 
* Sinkhorn optimal transport approach
* Multiscale version of sinkhorn transport
* Rpackage for sinkhorn transport (both multiscale and single scale, depends on gmra R package)
* Rpackage for optimal transport (both multiscale and single scale, requires any of CPLEX, GPLK or MOSEK, depends on gmra r package )

## Requirements

R packages:

* gmra https://github.com/suppechasper/gmra


C/C++ Libraries:

* Eigen
* CPLEX
     * Alternatively MOSEK or GLPK (currently setup for CPLEX but few modifications are required to change it to MOSEK or GLPK, edit the .cxx driver files in commandline accordingly to use the different solver header files and instantiate the solver you want to use. Then edit CMake files to point to the library specific header files and libraries. )

## R Package mop Installation


This package depends on linear programming libraries. Currently it supports GLPK,
CPLEX and MOSEK. CPLEX is typically the fastest.

To setup the package for a specific library edit src/Makevars and
src/mop_config.h in the mop directory

To install the package requires that all the libraries you are using are on R's linker path.
