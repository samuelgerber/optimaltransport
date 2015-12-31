This package depends on linear programming libraries. Currently it supports GLPK,
CPLEX and MOSEK. CPLEX is typically the fastest.

To setup the package for a specific library edit src/Makevars and
src/mop_config.h

To install the package requires that all the libraries you are using are on R's linker path.
