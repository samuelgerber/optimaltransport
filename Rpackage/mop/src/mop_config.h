//Which linear programming libraries should be made avaiable. Requires for R to
//have the corresponding libraries on its linker path. Currently all versions
//are active uncomment if you do not want to suport one of the libraries and
//adjust the linker arguments in Makevars.

//#define MOP_USE_CPLEX
//#define MOP_USE_MOSEK
#define MOP_USE_GLPK
