#ifndef MPS_H
#define MPS_H

#include "tsrproj.h"

using namespace std;

class MPS
{
    public:

        MPS (const int& _N, const int& _D, const bool& pbc, const int& CFLG=0);
        ~MPS ();

        inline TensorProj& getTsrProj (const int& i) {return mps[i];}
        inline const TensorProj& getTsrProj (const int& i) const {return mps[i];}

        void print ();

    private:

        int N;
        TensorProj* mps;

        //  initialize mps as an auxilliary function
        void initialize (const int& cfg);
};
#endif
