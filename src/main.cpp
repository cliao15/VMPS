#include "diag.h"
#include "vmps.h"

int main (int argc, char** argv)
{
    //  read Variable and Print
    RI ri;
    ri.readInputVar ();
    ri.printVal ();

    if (ri.MTH == 0)
    {
        //  using exact diagnalization method
        Diag qsDiagSolver;
        qsDiagSolver.MatrixH (ri);
        qsDiagSolver.DiagH ();
        qsDiagSolver.Output (ri);
    }
    else
    {
        //  using variational method
        VMPS qsVMPSsolver (ri);

        //  search for ground state
        cout << "\nVMPS solver for ground state -- begin:" << endl;
        qsVMPSsolver.minimizeE (0);
        cout << "\n\n\n\nVMPS solver for ground state -- end^^\n\n" << endl;

        //  search for first excited state
        if (ri.CALC_FES) 
        {
            cout << "\nVMPS solver for first excited state -- begin:" << endl;
            qsVMPSsolver.minimizeE (1);
            cout << "\n\n\n\nVMPS solver for first excited state -- end^^\n\n" << endl;
        }

        //  output res
        qsVMPSsolver.Output (ri);
    }

    return 0;
}
