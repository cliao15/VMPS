#ifndef DIAG_H
#define DIAG_H

#include <vector>
#include <ctime>
#include <iomanip>
#include <Eigen/Eigenvalues>

#include "readinput.h"

using namespace std;
using namespace Eigen;

#define TINY_Q 1E-10

class Diag
{
    public:
        
        Diag ();
        ~Diag ();

        //  calculate elements of Hamiltonian matrix
        void MatrixH (const RI& ri);

        //  calculate gs and/or fes properties
        void DiagH ();

        //  output
        int Output (const RI& ri);

        //  Energy and EigenVector
        VectorXd eigenval;
        MatrixXcd eigenvec;

    private:

        //  help calculate GS magnetization in X,Y,Z direction
        double Mxyz (const char&, const int&, const int&);

        //  help calculate GS correlation in X,Y,Z direaction
        double Cxyz (const char&, const int&, const int&, const int&, const bool&);

        //  cross pruduct
        void cross (const MatrixXcd&, const MatrixXcd&, MatrixXcd&) const;

        //  H elements due to n.n and n.n.n interaction
        void H_J (const char&, const int&, const int&, const double&, const bool&);

        //  H elements due to transverse and longitude field
        void H_M (const char&, const int&, const double&);

        //  Hamiltonian
        MatrixXcd H;

        //  Pauli Matrix
        MatrixXcd sx, sy, sz, id;

        //  time recorder
        double timeMH, timeDH, timeOUT;

};

#endif
