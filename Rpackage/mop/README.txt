This package depends on linear programming libraries. Currently it supports Lemon network simplex, GLPK,
CPLEX and MOSEK. Lemon and CPLEX network simplex are typically the fastest.

The default setup si lemon and is contained within the package. For using
another library edit src/Makevars and src/mop_config.h

To install the package requires that all the libraries you are using are on R's
linker path.
