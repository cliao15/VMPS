#include "mps.h"

MPS::MPS (const int& _N, const int& _D, const bool& pbc, const int& CFLG) 
    :N (_N), mps (NULL)
{
    mps = new TensorProj[N];
    if (pbc) {
        for (int i=0; i<N; i++) 
            mps[i].reAlloc (_D, _D);
    }
    else
    {
        if (N == 1) mps[0].reAlloc (1,1);
        else
        {
            mps[0].reAlloc (1, _D);
            mps[N-1].reAlloc (_D, 1);
            for (int i=1; i<N-1; i++)
                mps[i].reAlloc (_D, _D);
        }
    }


    initialize (CFLG); //initialization
}

MPS::~MPS () 
{
    if (mps != NULL)
    {
        delete [] mps;
        mps = NULL;
    }
}

void MPS::print ()
{
    for (int n=0; n<N; n++)
    {
        cout << "\nsite : " << n << "\n";
        for (int i=0; i<mps[n].getDtot (); i++)
        {
            cout << "(" << mps[n].getTensor (i).real ()
                 << ", " << mps[n].getTensor (i).imag ()
                 << ")\n";
        }
        cout << endl;
    }

    cout << "\n\n";
    return;
}

void MPS::initialize (const int& cfg)
{
    //  0: randomize configuration
    //  1: ferromagnetism
    //  2: anti-ferromagnetism

    if (cfg == 0)
    {
        //  randomization
        srand (time (NULL));    //  seed
        //  set tensor for each site
        for (int n=0; n<N; n++)
        {
            //  set DL*DR*d Matrix
            for (int i=0; i<mps[n].getDtot (); i++)
            {
                double real = (double)rand ()/RAND_MAX;
                mps[n].setTensor (i, complex<double>(real,0.0));
            }
        }
    }
    
    if (cfg == 1)
    {
        //  ferromagnetism
        for (int n=0; n<N; n++)
        {
            //  only the value at lefttop corner is one (d=1)
            //  otherwise zero 
            mps[n].setTensor (0,0,1,complex<double>(1.0,0.0));
        }
    }
    
    if (cfg == 2)
    {
        //  anti-ferromagnetism
        for (int n=0; n<N; n++)
        {
            if (n % 2 == 0) mps[n].setTensor (0,0,1,complex<double>(1.0,0.0));
            else mps[n].setTensor (0,0,0,complex<double>(1.0,0.0));
        }
    }

    return;
}



