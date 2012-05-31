#ifndef VMPS_H
#define VMPS_H

#include "mps.h"
#include "readinput.h"
#include <vector>
#include <map>
#include <string>
#include <complex>
#include <iomanip>

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

#define TINY 1E-10
#define USE_PREPARE 1

typedef map< string, complex<double> > CFGSTR;

class VMPS
{
    public:

        VMPS (const RI& ri);
        ~VMPS ();

        //  minimize Energy
        void minimizeE (const bool&);

        //  output
        int Output (const RI& ri); 

    private:

        //  set Hamiltonian
        void setMatrixH (const RI& ri);

        //  help calculate magnetization in x,y,z direaction
        double Mxyz (const char&, const bool&) const;

        //  help calculate correction (length r) in x,y,z direction
        double Cxyz (const char&, const int&, const bool&) const;

        //  storage 
        void initHstorage (const MPS&, MatrixXcd** Hstorage);
        void initCNstorage (const MPS&, const MPS&,  const MatrixXcd*, MatrixXcd*);
        
        //  projector in calculating excited states
        void calcprojector_onesite (const TensorProj&, const MatrixXcd&, const MatrixXcd&, MatrixXcd&);

        //  minimize E onesite
        double minimizeE_onesite (const int&, const bool&, const MatrixXcd&, TensorProj&);

        //  prepare site to be unitary tensor
        void prepare (MPS&, const bool&);
        void prepare_onesite (TensorProj&, MatrixXcd&, int&, const bool&);

        //  storage update
        void storage_update (const bool&, const int&, const TensorProj&, const bool&, const TensorProj&);

        //  output configuration
        void getConfig (const MPS&, CFGSTR&) const;

        int M, N; // M: total terms in Halmitonian, N: site
        MatrixXcd** hset;

        //  auxilliary storages 
        MatrixXcd** Hstorage;
        MatrixXcd*  Cstorage;
        MatrixXcd*  Nstorage;

        //  Pauli Matrix
        MatrixXcd sx, sy, sz, id;

        //  ground state / first exicited state
        double E0, E1;
        MPS* mpsGS;
        MPS* mpsFES;

        //  precision
        double tol;
        
        //  boundary condition
        bool pbc;

        //  time consumption
        double timeOUT, timeOPT, timeOPT_one;
};

#endif 

