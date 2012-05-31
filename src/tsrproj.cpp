#include "tsrproj.h"

TensorProj::TensorProj (const int& _D) :tensor (NULL), d(2)
{
    DL = DR = _D>1 ? _D:1;
    int Dtot = DL*DR*d;
    tensor = new complex<double> [Dtot];
    for (int i=0; i<Dtot; i++) tensor[i]=complex<double> (0.0,0.0);
}

TensorProj::~TensorProj () 
{
    if (tensor != NULL) 
    {
        delete [] tensor;
        tensor = NULL;
    }
}    

void TensorProj::reAlloc (const int& _DL, const int& _DR)
{
    //  free space 
    if (tensor != NULL) delete [] tensor;

    //  set Dims and Allocation
    DL = _DL; 
    DR = _DR;
    int Dtot = DL*DR*d;
    tensor = new complex<double> [Dtot];
    for (int i=0; i<Dtot; i++) tensor[i]=complex<double> (0.0,0.0);

    return;
}

TensorProj& TensorProj::operator = (const TensorProj& rhs)
{
    //  deallocate must be done first!
    if (tensor != NULL) 
    {
        delete [] tensor;
        tensor = NULL;
    }

    //  allocate
    DL = rhs.DL;
    DR = rhs.DR;

    int Dtot = DL*DR*d;
    tensor = new complex<double>[Dtot];

    //  copy data
    if (rhs.tensor != NULL)
    {
        for (int i=0; i<Dtot; i++) tensor[i]=rhs.tensor[i];
    }

    return *this;
}
