#ifndef READ_INPUT_H
#define READ_INPUT_H

#include "iostream"
#include "fstream"
#include "string"
#include "stdexcept"

using namespace std;

class RI 
{
    public:

        RI ():  NS(1), VD(1), MTH(0), INIT_CFG (0), CL(1),
                J1X(0.0), J1Y(0.0), J1Z(0.0), J2X(0.0), J2Y(0.0), J2Z(0.0),
                HX(0.0), HY(0.0), HZ(0.0), TOL(5e-3), CALC_FES(false), 
                OUT_CFG (false), USE_PBC(false), OUT_MAG (0), OUT_CORR (0)
        {}
       ~RI () {}

        //  read variables from input file
        int readInputVar ();

        //  print varaibles read from input file
        int printVal ();

        //  variables 
        int NS, VD, MTH, INIT_CFG, CL;
        double J1X, J1Y, J1Z, J2X, J2Y, J2Z, HX, HY, HZ, TOL;
        bool CALC_FES, OUT_MAG, OUT_CORR, OUT_CFG, USE_PBC;
};

#endif
