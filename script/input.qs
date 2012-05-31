#   ------------------------------------------------------------------------------------------------------------------
#   input file to program qsANNNI - simulator of quantum state of 1D ANNNI model
#   reasonable values are strongly recommended otherwise unexpected errors may occur!!
#
#   Two methods are available 
#   (1) exact diagnalization (give all eigenstates)
#   (2) variational Matrix Product State (MPS) (Ground State (GS) is default, First Exicted State (FES) is optional)
#
#   Hamiltonian = -sum{j=1,N}
#   (
#       J1X*sx(j)*sx(j+1)+J1Y*sy(j)*sy(j+1)+J1Z*sz(j)*sz(j+1)
#       +
#       J2X*sx(j)*sx(j+2)+J2Y*sy(j)*sy(j+2)+J2Z*sz(j)*sz(j+2)
#       +
#       HX*sx(j)+HY*sy(j)+HZ*sz(j)
#       )
#   both Open Boundary Condition (OBC) or Periodical Boundary Condition (PBC) is optional
#
#   Output options:
#   (1) Magnetization in X,Y,Z direction
#   (2) Correlation in X,Y,Z direction (set correlation length CL, e.g. sx(j)*sx(j+CL))
#   (3) Configurations
#   -----------------------------------------------------------------------------------------------------------------

#   number of sites
#   default :1
NS      7

#   dimension of virtual bond
#   default :1
VD      1

#   simulation method:
#   Exact Diagonalization :0
#   Variational MPS :1
#   default :0
MTH     0

#   initial configuration (VMPS only)
#   random cfg :0
#   ferromagnetism :1
#   anti-ferromagnetism :2
#   default :0
INIT_CFG 0

#   length of corrlation (used for output Cxx, Cyy, Czz)
#   default :1
CL      1

#   nearest neighbor interaction strength in X direction
#   default :0.0
J1X      1

#   nearest neighbor interaction strength in Y direction
#   default :0.0
J1Y      0

#   nearest neighbor interaction strength in Z direction
#   default :0.0
J1Z     0.0

#   next nearest neighbor interaction strength in X direction
#   default :0.0
J2X     0.0

#   next nearest neighbor interaction strength in Y direction
#   default :0.0
J2Y     0.0

#   next nearest neighbor interaction strength in Z direction
#   default :0.0
J2Z     0.0

#   field strength in X direction
#   default :0.0
HX      1e-3

#   field strength in Y direction
#   default :0.0
HY      0.0

#   field strength in Z direction
#   default :0.0
HZ      0.02

#   error tolerance
#   default :5e-3
TOL     1e-5 

#   calculate first exicited state? (VMPS only)
#   default :0
CALC_FES     0

#   use periodic boundary condition otherwise open boundary condition will be used
#   default :0
USE_PBC      1

#   output Magnetization in (X,Y,Z) direction
#   default :0
OUT_MAG      1

#   output Correlation in (X,Y,Z) direaction
#   default :0
OUT_CORR     0

#   output configurations of ground state
#   default :0
OUT_CFG      0
