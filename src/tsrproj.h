#ifndef TENSOR_PROJECTOR_H
#define TENSOR_PROJECTOR_H

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <complex>
#include <cassert>

using namespace std;

class TensorProj
{
    public:

        //  constructor
        TensorProj (const int& _D=1);
        
        //  deconstructor
        ~TensorProj ();

        //  overload = operator
        TensorProj& operator = (const TensorProj&);

        /**
         * get dimension of both physical system as well as Auxilliary systems
         */
        inline int getDL () const {return DL;}
        inline int getDR () const {return DR;}
        inline int getd  () const {return d; }
        inline int getDtot () const {return DL*DR*d;}

        //  resize and allocation
        //  alloc tensor based on virtual and physical indice
        void reAlloc (const int&, const int&);

        //  read
        inline complex <double> getTensor (const int& i) const {return tensor[i];}
        inline complex <double> getTensor (
                const int& _dl, const int& _dr, const int& _d
                ) const  {
            return tensor[getIndex (_dl,_dr,_d)];
        } 

        //  write
        inline void setTensor (
                const int& _dl, const int& _dr, const int& _d, const complex<double>& val
                ) {
            tensor[getIndex (_dl,_dr,_d)] = val;
        }

        inline void setTensor (const int& index, const complex<double>& val) {
            assert (index >= 0); tensor[index] = val;
        }

    private:

        //  get combined indice from (l,r,k)
        inline int getIndex (
                const int& _dl, const int& _dr, const int& _d
                ) const {
            return  _dl*DR*d + _dr*d + _d;
        }

        //  dim of virtual bonds
        //  dim of physics index is fixed for quantum spin system
        int DL, DR, d;

        /**
         * compress 3 indices into 1-dim array
         */
        complex <double> * tensor;
};

#endif 
