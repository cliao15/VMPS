#include "vmps.h"

VMPS::VMPS (const RI& ri)
    : E0 (0.0), E1 (0.0)
      , Hstorage(NULL), Cstorage(NULL), Nstorage(NULL)
      , mpsGS (NULL), mpsFES (NULL)
      , timeOUT (0.0), timeOPT (0.0), timeOPT_one (0.0)
{
    //  auxilliary matrix
    sx = sy = sz = id = MatrixXcd::Zero (2,2);
    sx (1,0) = sx (0,1) = 1.0;
    sy (0,1) = complex<double>(0.0,-1.0); sy (1,0) = complex<double>(0.0,1.0);
    sz (0,0) = 1.0; sz (1,1) = -1.0;
    id (0,0) = id (1,1) = 1.0;
    tol = ri.TOL;
    pbc = ri.USE_PBC;

    //  hset
    setMatrixH (ri);

    //  init MPS
    mpsGS = new MPS (ri.NS, ri.VD, ri.USE_PBC, ri.INIT_CFG);
    if (ri.CALC_FES) mpsFES = new MPS (ri.NS, ri.VD, ri.USE_PBC, ri.INIT_CFG);
}

VMPS::~VMPS ()
{
    if (Hstorage != NULL)
    {
        for (int i=0; i<M; i++) delete [] Hstorage[i];
        delete [] Hstorage;
        Hstorage = NULL;
    }

    if (Cstorage != NULL)
    {
        delete [] Cstorage;
        Cstorage = NULL;
    }

    if (Nstorage != NULL)
    {
        delete [] Nstorage;
        Nstorage = NULL;
    }

    if (hset != NULL)
    {
        for (int i=0; i<M; i++) delete [] hset[i];
        delete [] hset;
        hset = NULL;
    }

    if (mpsGS != NULL)
    {
        delete mpsGS;
        mpsGS = NULL;
    }

    if (mpsFES != NULL)
    {
        delete mpsFES;
        mpsFES = NULL;
    }
}

void VMPS::setMatrixH (const RI& ri)
{
    //  copy data
    N = ri.NS;

    /****************
     * TEST BEGIN
     * *************/

    /**
    M = 2*N-1;
    hset = new MatrixXcd*[M];
    for (int m=0; m<M; m++)
    {
        hset[m] = new MatrixXcd[N];
        for (int j=0; j<N; j++)
        {
            hset[m][j] = id;
        }
    }

    for (int j=0; j<N-1; j++)
    {
        hset[2*j][j] = -sx;
        hset[2*j][j+1] = sx;
        hset[2*j+1][j] = -sz;
    }

    hset[M-1][N-1] = -sz;
    return;
    */

    /****************
     * TEST END
     * *************/


    //  calculate M
    M = 0;
    if (fabs (ri.J1X) > TINY) M += (ri.USE_PBC==true) ? N:N-1;
    if (fabs (ri.J1Y) > TINY) M += (ri.USE_PBC==true) ? N:N-1;
    if (fabs (ri.J1Z) > TINY) M += (ri.USE_PBC==true) ? N:N-1;
    if (fabs (ri.J2X) > TINY) M += (ri.USE_PBC==true) ? N:N-2;
    if (fabs (ri.J2Y) > TINY) M += (ri.USE_PBC==true) ? N:N-2;
    if (fabs (ri.J2Z) > TINY) M += (ri.USE_PBC==true) ? N:N-2;
    if (fabs (ri.HX)  > TINY) M += N;
    if (fabs (ri.HY)  > TINY) M += N;
    if (fabs (ri.HZ)  > TINY) M += N;

    //  allocate space
    hset = new MatrixXcd*[M];
    for (int i=0; i<M; i++) 
        hset[i] = new MatrixXcd[N];

    //  assignment
    int p = 0;

    //  J1X
    if (fabs (ri.J1X) > TINY) 
    {
        for (int i=0; i<N-1; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J1X*sx;
                else if (j==i+1) hset[p+i][j] = sx;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-1][0] = -ri.J1X*sx;
            hset[p+N-1][N-1] = sx;
            for (int j=1; j<N-1; j++) hset[p+N-1][j] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-1;
    }

    //  J1Y
    if (fabs (ri.J1Y) > TINY) 
    {
        for (int i=0; i<N-1; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J1Y*sy;
                else if (j==i+1) hset[p+i][j] = sy;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-1][0] = -ri.J1Y*sy;
            hset[p+N-1][N-1] = sy;
            for (int j=1; j<N-1; j++) hset[p+N-1][j] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-1;
    }

    //  J1Z
    if (fabs (ri.J1Z) > TINY) 
    {
        for (int i=0; i<N-1; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J1Z*sz;
                else if (j==i+1) hset[p+i][j] = sz;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-1][0] = -ri.J1Z*sz;
            hset[p+N-1][N-1] = sz;
            for (int j=1; j<N-1; j++) hset[p+N-1][j] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-1;
    }

    //  J2X
    if (fabs (ri.J2X) > TINY)
    {
        for (int i=0; i<N-2; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J2X*sx;
                else if (j==i+2) hset[p+i][j] = sx;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-2][0] = -ri.J2X*sx;
            hset[p+N-2][N-2] = sx;
            for (int j=1; j<N-2; j++) hset[p+N-2][j] = id;
            hset[p+N-2][N-1] = id;

            hset[p+N-2][1] = -ri.J2X*sx;
            hset[p+N-2][N-1] = sx;
            for (int j=2; j<N-1; j++) hset[p+N-1][j] = id;
            hset[p+N-1][0] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-2;
    }

    //  J2Y
    if (fabs (ri.J2Y) > TINY)
    {
        for (int i=0; i<N-2; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J2Y*sy;
                else if (j==i+2) hset[p+i][j] = sy;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-2][0] = -ri.J2Y*sy;
            hset[p+N-2][N-2] = sy;
            for (int j=1; j<N-2; j++) hset[p+N-2][j] = id;
            hset[p+N-2][N-1] = id;

            hset[p+N-2][1] = -ri.J2Y*sy;
            hset[p+N-2][N-1] = sy;
            for (int j=2; j<N-1; j++) hset[p+N-1][j] = id;
            hset[p+N-1][0] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-2;
    }

    //  J2Z
    if (fabs (ri.J2Z) > TINY)
    {
        for (int i=0; i<N-2; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j==i) hset[p+i][j] = -ri.J2Z*sz;
                else if (j==i+2) hset[p+i][j] = sz;
                else hset[p+i][j] = id;
            }
        }

        if (ri.USE_PBC)
        {
            hset[p+N-2][0] = -ri.J2Z*sz;
            hset[p+N-2][N-2] = sz;
            for (int j=1; j<N-2; j++) hset[p+N-2][j] = id;
            hset[p+N-2][N-1] = id;

            hset[p+N-2][1] = -ri.J2Z*sz;
            hset[p+N-2][N-1] = sz;
            for (int j=2; j<N-1; j++) hset[p+N-1][j] = id;
            hset[p+N-1][0] = id;
        }

        //  update
        p += (ri.USE_PBC==true) ? N:N-2;
    }

    //  Hx
    if (fabs (ri.HX) > TINY)
    {
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j == i) hset[p+i][j] = -ri.HX*sx;
                else hset[p+i][j] = id;
            }
        }

        //  update
        p += N;
    }

    //  Hy
    if (fabs (ri.HY) > TINY)
    {
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j == i) hset[p+i][j] = -ri.HY*sy;
                else hset[p+i][j] = id;
            }
        }

        //  update
        p += N;
    }

    //  Hz
    if (fabs (ri.HZ) > TINY)
    {
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<N; j++)
            {
                if (j == i) hset[p+i][j] = -ri.HZ*sz;
                else hset[p+i][j] = id;
            }
        }

        //  update
        p += N;
    }

    //  fatal assertion
//    cout << "\nM = " << M << endl;
    assert (p == M);

    return;
}

void VMPS::initHstorage (
        const MPS& mps, 
        MatrixXcd** Hstorage // dim M*(N+1)
        )
{
    //  for site 0 and site N
    for (int m=0; m<M; m++) 
    {
        //  n=0;
        int dim0 = mps.getTsrProj(0).getDL ();
   //     cout << "\ndim0 = " << dim0 << endl;
        dim0 *= dim0;
        Hstorage[m][0] = MatrixXcd::Identity (dim0, dim0);

        //  n=N;
        int dimN = mps.getTsrProj(N-1).getDR ();
        dimN *= dimN;
        Hstorage[m][N] = MatrixXcd::Identity (dimN, dimN);
    }
    
    for (int n=N-1; n>0; n--)
    {
        for (int m=0; m<M; m++)
        {
            /**
             * (1) contract h[m][n] wity proj A and proj A+
             */
            const TensorProj& proj = mps.getTsrProj (n);

            const int DL = proj.getDL ();
            const int DR = proj.getDR ();
            const int d  = proj.getd  ();

            MatrixXcd A = MatrixXcd::Zero (DL*DR,d);
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        A (i*DR+j,k) = proj.getTensor (i,j,k);
                    }
                }
            }

            assert (hset[m][n].rows () == d);
            MatrixXcd B = A * hset[m][n] *A.adjoint (); // contraction

            /**
             * (2) reshape B to matrix C with dim (DL*DL, DR*DR)
             */
            MatrixXcd C = MatrixXcd::Zero (DL*DL, DR*DR);
            for (int i1=0; i1<DL; i1++)
            {
                for (int i2=0; i2<DR; i2++)
                {
                    for (int j1=0; j1<DL; j1++)
                    {
                        for (int j2=0; j2<DR; j2++)
                        {
                            C(i1*DL+j1, i2*DR+j2) = B (i1*DR+i2, j1*DR+j2);
                        }
                    }
                }
            }
            
    //        if (m == 2 && n == 1) cout << "\nC = " << C << endl;
            
            /**
             * (3) contract C with Hstorage (n+1)
             */
            Hstorage[m][n] = C*Hstorage[m][n+1];
            /*
            cout << "\nm = " << m << " n = " << n 
                 << " hset = \n" << hset[m][n]  << endl 
                 << " Hstorage = \n" << Hstorage[m][n] << endl;
                */
        }
    }

    //  end func initHstorage
    return;
}

void VMPS::initCNstorage (
        const MPS& mpsA,     //  state needed to be solved 
        const MPS& mpsB,     //  referenced state
        const MatrixXcd* X,  //  operator
        MatrixXcd* CNstorage  //  N+1 element!
        )
{

    //  site 0
    int dim0 = mpsA.getTsrProj(0).getDL() * mpsB.getTsrProj(0).getDL();
    CNstorage[0] = MatrixXcd::Identity (dim0, dim0);

    //  site N
    int dimN = mpsA.getTsrProj(N-1).getDR() * mpsB.getTsrProj(N-1).getDR();
    CNstorage[N] = MatrixXcd::Identity (dimN, dimN);
    
    for (int n=N-1; n>0; n--)
    {
        /**
         * (1) contract A with X (n) and B+ <B|X|A>
         */

        //  A
        const TensorProj& projA = mpsA.getTsrProj (n);

        const int DLA = projA.getDL ();
        const int DRA = projA.getDR ();
        const int dA  = projA.getd  ();

        MatrixXcd A = MatrixXcd::Zero (DLA*DRA,dA);
        for (int i=0; i<DLA; i++)
        {
            for (int j=0; j<DRA; j++)
            {
                for (int k=0; k<dA; k++)
                {
                    A (i*DRA+j,k) = projA.getTensor (i,j,k);
                }
            }
        }

        //  B
        const TensorProj& projB = mpsB.getTsrProj (n);

        const int DLB = projB.getDL ();
        const int DRB = projB.getDR ();
        const int dB  = projB.getd  ();
        
        MatrixXcd B = MatrixXcd::Zero (dB, DLB*DRB);
        for (int i=0; i<DLB; i++)
        {
            for (int j=0; j<DRB; j++)
            {
                for (int k=0; k<dB; k++)
                {
                    B (k, i*DRB+j) = projB.getTensor (i,j,k);
                }
            }
        }

        MatrixXcd C (DLA*DRA, DLB*DRB);
        if (X !=NULL) C = A*X[n]*B.conjugate ();
        else C = A*B.conjugate ();

        /**
         * (2) reshape C to be a matrix with dim (DLA*DLB, DRA*DRB)
         */
        MatrixXcd D = MatrixXcd::Zero (DLA*DLB, DRA*DRB);
        for (int i1=0; i1<DLA; i1++)
        {
            for (int i2=0; i2<DRA; i2++)
            {
                for (int j1=0; j1<DLB; j1++)
                {
                    for (int j2=0; j2<DRB; j2++)
                    {
                        D (i1*DLB+j1, i2*DRB+j2) = C (i1*DRA+i2, j1*DRB+j2);
                    }
                }
            }
        }

        /**
         * (3) contract C with Cstorage (n+1)
         */
        CNstorage[n] = D*CNstorage[n+1];

//        cout << "\nn = " << n << " CNstorage = \n" << CNstorage[n] << endl;
    }

    //  end func initCstorage
    return;
}

void VMPS::calcprojector_onesite (
        const TensorProj& A,
        const MatrixXcd& Cleft,
        const MatrixXcd& Cright,
        MatrixXcd& Proj
        )
{
    //  A
    const int DL = A.getDL ();
    const int DR = A.getDR ();
    const int d  = A.getd  ();

    /**
     * trace Cleft with Cright
     */
    MatrixXcd B = Cleft.transpose () * Cright.transpose ();

    /**
     * reshape B with dim (yL*yR, DL*DR)
     */
    const int yL = B.rows ()/DL;
    const int yR = B.cols ()/DR;
    const int dy = d;

    MatrixXcd C = MatrixXcd::Zero (yL*yR, DL*DR);
    for (int i=0; i<yL; i++)
    {
        for (int j=0; j<DL; j++)
        {
            for (int k=0; k<yR; k++)
            {
                for (int l=0; l<DR; l++)
                {
                    C (i*yR+k, j*DR+l) = B (i*DL+j, k*DR+l);
                }
            }
        }
    }

    /**
     * reshape A with dim (DL*DR, d)
     */
    MatrixXcd D = MatrixXcd::Zero (DL*DR,d);
    for (int i=0; i<DL; i++)
    {
        for (int j=0; j<DR; j++)
        {
            for (int k=0; k<d; k++)
            {
                D (i*DR+j,k) = A.getTensor (i,j,k);
            }
        }
    }

    /**
     * contract C with A+
     */
    MatrixXcd E = C * D.conjugate (); 

    //  reshape E into 1D vector Vec
    MatrixXcd Vec = MatrixXcd::Zero (yL*yR*dy,1);
    for (int i=0; i<yL; i++)
    {
        for (int j=0;j<yR; j++)
        {
            for (int k=0; k<dy; k++)
            {
                Vec (i*yR*dy+j*dy+k,0) = E (i*yR+j,k);
            }
        }
    }

    /**
     * find projectors
     */
    const int sizey = yL*yR*dy;
    MatrixXcd orth (sizey, sizey+1);
    for (int i=0; i<sizey; i++) 
    {
        for (int j=0; j<sizey+1; j++)
        {
            if (j == 0) orth (i,0) = Vec (i,0);
            else if (j == i+1) orth (i,i+1) = 1.0;
            else orth (i,j) = 0.0;
        }
    }

    //  schdmit orthgonalization
    VectorXcd vxd = VectorXcd::Zero (sizey); 
    VectorXcd w = VectorXcd::Zero (sizey);

    vector<int> count;  //count the non-zero basis

    for (int j=0; j<=sizey; j++)
    {
        //  fetch data
        for (int i=0; i<sizey; i++) vxd (i) = orth (i,j);

        //  schdmit process
        for (int k=j-1; k>=0; k--)
        {
            for (int i=0; i<sizey; i++) w (i) = orth (i,k);
            double wdot = w.dot(w).real();
            if (wdot > TINY) vxd = vxd - vxd.dot(w)/w.dot(w)*w;
        }
        
        //  normalize
        double norm = sqrt (vxd.dot(vxd).real ());
        if (norm > TINY) vxd = vxd / norm;

        //  update
        for (int i=0; i<sizey; i++) orth (i,j) = vxd (i);

        //  write to proj
        //  ATTENTION:  ZERO vectors must be omited
        if (j > 0 && sqrt(vxd.dot(vxd).real()) > TINY) count.push_back (j);
    }

    //  copy data to proj
    Proj.resize (sizey, count.size());
    for (int j=0; j<count.size (); j++)
    {
        for (int i=0; i<sizey; i++)
        {
            Proj (i,j) = orth (i,count[j]);
        }
    }

//  cout << "\ncount = " << count.size() << " Proj = \n" << Proj << endl;
    
    //  end calcproj_onesite
    return;
}

double VMPS::minimizeE_onesite (
        const int& sj,   //  site j
        const bool& gs, //  0 for ground state and 1 for first excited state
        const MatrixXcd& Proj,    //  projector
        TensorProj& A    //  optimized tensor
        )
{
    //  TIME
    time_t start, end;
    time (&start);

    /**
     * calculate Heff and Neff
     */

    //  ATTENTION! DAl may not be equal to A.getDL () since
    //  we have changed DAl by performing SVD
    const int DAl = sqrt (Hstorage[0][sj].cols ());
    const int DAr = sqrt (Hstorage[0][sj+1].rows ());
    const int d   = A.getd  ();
    const int sizeA = DAl * DAr * d;

    //  calculate Heff
    MatrixXcd Heff = MatrixXcd::Zero (sizeA, sizeA);
    for (int m=0; m<M; m++)
    {
        //  trace Hstorage[m][sj] and trace Hstorage[m][sj+1]
        MatrixXcd tmpTrH = Hstorage[m][sj].transpose () * Hstorage[m][sj+1].transpose ();

        //  temporary H storage for each m
        MatrixXcd tmp = MatrixXcd::Zero (sizeA, sizeA);

        //  assignment
        for (int i1=0; i1<DAl; i1++)
        {
            for (int i2=0; i2<DAl; i2++)
            {
                for (int j1=0; j1<DAr; j1++)
                {
                    for (int j2=0; j2<DAr; j2++)
                    {
                        for (int k1=0; k1<d; k1++)
                        {
                            for (int k2=0; k2<d; k2++)
                            {
                                tmp (i1*DAr*d+j1*d+k1, i2*DAr*d+j2*d+k2) = 
                                    tmpTrH (i1*DAl+i2, j1*DAr+j2) * hset[m][sj](k1,k2);
                            }
                        }
                    }
                }
            }
        }

        //  addition
        Heff += tmp;
    }

    //  Neff is needed for PBC
    MatrixXcd Neff = MatrixXcd::Zero (sizeA, sizeA);

    if (pbc || !USE_PREPARE)
    {
        //  calculate Neff

        //  trace Hstorage[m][sj] and trace Hstorage[m][sj+1]
        MatrixXcd tmpTrN = Nstorage[sj].transpose () * Nstorage[sj+1].transpose ();

        for (int i1=0; i1<DAl; i1++)
        {
            for (int i2=0; i2<DAl; i2++)
            {
                for (int j1=0; j1<DAr; j1++)
                {
                    for (int j2=0; j2<DAr; j2++)
                    {
                        for (int k=0; k<d; k++)
                        {
                            Neff (i1*DAr*d+j1*d+k, i2*DAr*d+j2*d+k) = tmpTrN (i1*DAl+i2, j1*DAr+j2);
                        }
                    }
                }
            }
        }
    }
    
    //  projects on orthognal subspace
    if (gs) 
    {
        Heff = Proj.adjoint () * Heff * Proj;
        if (pbc) Neff = Proj.adjoint () * Neff * Proj;
    }
    
    Heff = (Heff + Heff.adjoint ())/2.0; // enforce hermiticity
    if (pbc || !USE_PREPARE) Neff = (Neff + Neff.adjoint ())/2.0; // enforce hermiticity

    double lambda = 0.0;
    VectorXcd v = VectorXcd::Identity (sizeA,1);

    if (pbc || !USE_PREPARE)
    {
        /**
         * improve the condition of Neff
         */
    //    cout << "\nNeff.row = " << Neff.rows () << " col = " << Neff.cols () << endl;

        SelfAdjointEigenSolver<MatrixXcd> es1(Neff);
        int eignum = es1.eigenvectors ().cols ();

        //  find number of non-zero eigenvalues;
        vector<int> ind;
        for (int i=0; i<eignum; i++) 
            if (fabs(es1.eigenvalues ()[i]) > TINY) 
                ind.push_back (i);
    //    cout << "\nnumber of non-zero eigenvalues = " << ind.size () << endl;

        if (ind.size () < eignum)
        {
            //  find projector Proj2
            MatrixXcd Proj2 (eignum, ind.size());
            for (int i=0; i<ind.size(); i++)
                Proj2.col(i) = es1.eigenvectors ().col(ind[i]);

            //  obtain projected Heff and Neff
            MatrixXcd Heffp = Proj2.adjoint () * Heff * Proj2;
            MatrixXcd Neffp = Proj2.adjoint () * Neff * Proj2;

            //  solver for generalized eigenvalue
            GeneralizedSelfAdjointEigenSolver<MatrixXcd> es2(Heffp, Neffp);

            //  eigenvalue
            lambda = es2.eigenvalues ()[0];
            //         cout << "\neigenvalues are :\n" << es.eigenvalues () << endl;

            //  eigenvector
            if (gs) v = Proj * Proj2 * es2.eigenvectors().col(0);
            else v = Proj2 * es2.eigenvectors().col(0);
        }
        else
        {
            //  solver for generalized eigenvalue
            GeneralizedSelfAdjointEigenSolver<MatrixXcd> es2(Heff, Neff);

            //  eigenvalue
            lambda = es2.eigenvalues ()[0];
            //         cout << "\neigenvalues are :\n" << es.eigenvalues () << endl;

            //  eigenvector
            if (gs) v = Proj * es2.eigenvectors().col(0);
            else v = es2.eigenvectors().col(0);
        }
    }
    else
    {
        //  solver for eigenvalue
        SelfAdjointEigenSolver<MatrixXcd> eigensolver(Heff);

        //  eigenvalue
        lambda = eigensolver.eigenvalues () (0,0);

        //  eigenvector
        if (gs) v = Proj * eigensolver.eigenvectors ().col(0);
        else v = eigensolver.eigenvectors ().col (0);
    }

    //  nan 
    if (lambda != lambda)
    {
        cerr << "\n\nBad Conditional Number -- Nan FOUND!\n\n" << endl;
        exit (1);
    }

 //   cout << "\nlambda = " << lambda << " | v = \n" << v << endl;

    /**
     * write to tensor
     */
    A.reAlloc (DAl, DAr);

    for (int i=0; i<DAl; i++)
    {
        for (int j=0; j<DAr; j++)
        {
            for (int k=0; k<d; k++)
            {
                A.setTensor (i, j, k, v(i*DAr*d+j*d+k));
            }
        }
    }

    //  time
    time (&end);
    timeOPT_one += difftime (end, start); 

    return lambda;
    //  end
}


void VMPS::minimizeE (const bool& gs) //  0 for ground state, 1 for excited state
{
    //  time
    time_t start, end;
    time (&start);

    //  target MPS
    MPS* mpsTgt = NULL;
    if (gs) mpsTgt = mpsFES;
    else mpsTgt = mpsGS;

    /**
     * prepare a mps using gauge transformation such that Neff is the identity for the first spin
     */
    if (USE_PREPARE) prepare (*mpsTgt, true);

    /**
     * storage initialization
     */

    //  Hstorage
    Hstorage = new MatrixXcd*[M];
    for (int i=0; i<M; i++) Hstorage[i] = new MatrixXcd[N+1];
    initHstorage (*mpsTgt, Hstorage);

    if (pbc || !USE_PREPARE)
    {
        //  Nstorage
        Nstorage = new MatrixXcd[N+1];
        initCNstorage (*mpsTgt, *mpsTgt, NULL, Nstorage);
    }

    //  Cstorage
    if (gs)
    {
        Cstorage = new MatrixXcd[N+1];
        initCNstorage (*mpsFES, *mpsGS, NULL, Cstorage);
    }

    /**
     * VMPS method by sweeping and optimizing each tensor
     */
    int count = 0;
    while (1)
    {
        cout << "\n\nround :" << count << endl;

        vector<double> Evalues;

        /************************************
         * cycle 1: j -> j+1 (from 1 to N-1)
         * **********************************/
        cout << "\nfrom left to right :";
        for (int j=0; j<N-1; j++)
        {
            //  projector calculation
            MatrixXcd P;    //  P will be resized in function calcprojector_onesite
            if (gs)  calcprojector_onesite ((*mpsGS).getTsrProj(j), Cstorage[j], Cstorage[j+1], P);

            //  optimization
            double E = minimizeE_onesite (j, gs, P, (*mpsTgt).getTsrProj (j));
            cout << "\nsite : " << j << " | energy = " << setiosflags (ios::fixed) << setprecision (10) << E;

            //  prepare_onesite and update mps{j}
            int DB = -1;
            MatrixXcd S;
            if (USE_PREPARE) prepare_onesite ((*mpsTgt).getTsrProj (j), S, DB, false); 

            //  write Evalues
            Evalues.push_back (E);

            //  storage update
            if (gs) storage_update (false, j, (*mpsFES).getTsrProj (j), true, (*mpsGS).getTsrProj (j));
            else storage_update (false, j, (*mpsGS).getTsrProj (j), false, (*mpsGS).getTsrProj (j));
        }

	    /************************************
         * cycle 1: j -> j-1 (from N to 2)
         * **********************************/

        //  since it is possible the loop will break after converging to ground state,
        //  we have to store U for the first tensor 
        MatrixXcd U;

        cout << "\n\nfrom right to left :";
        for (int j=N-1; j>0; j--)
        {
            //  projector calculation
            MatrixXcd P;    //  P will be resized in function calcprojector_onesite
            if (gs)  calcprojector_onesite ((*mpsGS).getTsrProj(j), Cstorage[j], Cstorage[j+1], P);

            //  optimization
            double E = minimizeE_onesite (j, gs, P, (*mpsTgt).getTsrProj (j));
            cout << "\nsite : " << j << " | energy = " << E;

            //  prepare_onesite and update mps{j}
            int DB = -1;
            if (USE_PREPARE) prepare_onesite ((*mpsTgt).getTsrProj (j), U, DB, true); 

            //  write Evalues
            Evalues.push_back (E);

            //  storage update
            if (gs) storage_update (true, j, (*mpsFES).getTsrProj (j), true, (*mpsGS).getTsrProj (j));
            else storage_update (true, j, (*mpsGS).getTsrProj (j), false, (*mpsGS).getTsrProj (j));
        }

        /**
         * evaluate the qualification of current tensor
         */
        
        // calculate expectation
        double mean = 0;
        for (int i=0; i<Evalues.size(); i++) mean += Evalues[i];
        mean /= Evalues.size ();

        // standard deviation
        double std = 0;
        for (int i=0; i<Evalues.size(); i++) std += (Evalues[i]-mean)*(Evalues[i]-mean);
        std = sqrt (std/Evalues.size());
        
        /**
         * if converge, then exit
         */
        if (std/fabs(mean) < tol)
        {
            if (USE_PREPARE)
            {
                TensorProj& tsr0 = (*mpsTgt).getTsrProj (0);

                //  contract mps[0] with U
                const int DL0 = tsr0.getDL ();
                const int DR0 = tsr0.getDR ();
                const int d0  = tsr0.getd  ();

                MatrixXcd tmp1 = MatrixXcd::Zero (DL0*d0, DR0);
                for (int i=0; i<DL0; i++)
                {
                    for (int j=0; j<DR0; j++)
                    {
                        for (int k=0; k<d0; k++)
                        {
                            tmp1 (i*d0+k, j) = tsr0.getTensor (i,j,k);
                        }
                    }
                }

                MatrixXcd tmp2 = tmp1 * U;

                //  rewrite to projector 0
                tsr0.reAlloc (DL0, U.cols ());
                for (int i=0; i<DL0; i++)
                {
                    for (int j=0; j<U.cols(); j++)
                    {
                        for (int k=0; k<d0; k++)
                        {
                            tsr0.setTensor (i, j, k, tmp2 (i*d0+k,j));
                        }
                    }
                }
            }

            //  ground state E0
            //  first excited state E1
            if (gs) E1 = Evalues.back ();
            else    E0 = Evalues.back ();

            break;
        }

        count ++;
    }

    /**
     * free space Hstorage, Nstorage, Cstorage
     */
    if (Hstorage != NULL)
    {
        //  Hstorage M*N
        for (int i=0; i<M; i++) delete [] Hstorage[i];
        delete [] Hstorage;
        Hstorage = NULL;
    }

    if (Nstorage != NULL)
    {
        //  Nstorage N
        delete [] Nstorage;
        Nstorage = NULL;
    }
    if (Cstorage != NULL)
    {
        //  Cstorage N
        delete [] Cstorage;
        Cstorage = NULL;
    }

    //  time
    time (&end);
    timeOPT += difftime (end, start);

    //  end VMPS
    return;
}

void VMPS::prepare  (MPS& mps, const bool& dir)
{
    if (dir == false)
    {
        //  prepare from left to right
        for (int n=0; n<N-1; n++)   //  omit the last site
        {
            int DB = -1;
            MatrixXcd U;
            prepare_onesite (mps.getTsrProj (n),U,DB,false);

            //  U: DB * DR
            //  fold tensor to matrix T of dimension: DL * (DR*d) 
            //  and then perform U*T and then unfold the result to tensor DB * DR * d
            TensorProj& tsr = mps.getTsrProj (n+1);

            const int DL = tsr.getDL ();
            const int DR = tsr.getDR ();
            const int d  = tsr.getd  ();

            MatrixXcd T = MatrixXcd::Zero (DL, DR*d);
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        T (i,j*d+k) = tsr.getTensor (i,j,k);
                    }
                }
            }

            //  contraction
            MatrixXcd W = U*T;    //  dim: DB * (d*DR)

            //  rewrite to tensor project
            tsr.reAlloc (DB, DR);
            for (int i=0; i<DB; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        tsr.setTensor (i, j, k, W (i, j*d+k));
                    }
                }
            }
        }

        return;
    }
    else
    {
        //  prepare from right to left
        for (int n=N-1; n>0; n--)   //  omit the first site
        {
            int DB = -1;
            MatrixXcd U;
            prepare_onesite (mps.getTsrProj (n), U, DB, true);

            //  U: DL * DB 
            //  fold tensor to matrix T of dimension: (DL*d) * DR
            //  and then perform T*U and then unfold the result to tensor DL * DB * d
            TensorProj& tsr = mps.getTsrProj (n-1);

            const int DL = tsr.getDL ();
            const int DR = tsr.getDR ();
            const int d  = tsr.getd  ();

            MatrixXcd T = MatrixXcd::Zero (DL*d, DR);
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        T (i*d+k,j) = tsr.getTensor (i,j,k);
                    }
                }
            }

            //  contraction
            MatrixXcd W = T*U;    //  dim: (DL*d) * DB

            //  new tensor
            tsr.reAlloc (DL, DB);
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DB; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        tsr.setTensor (i, j, k, W (i*d+k,j));
                    }
                }
            }
        }

        return;
    }

    //  end
    return;
}

void VMPS::prepare_onesite (TensorProj& tsrA, MatrixXcd& S, int& DB, const bool& dir)
{
    const int DL = tsrA.getDL ();
    const int DR = tsrA.getDR ();
    const int d  = tsrA.getd  ();

    if (dir == false)
    {
        //  prepare from left to right

        /**
         * construct matrix fold with dim (d*DL, DR)
         */
        MatrixXcd fold = MatrixXcd::Zero (d*DL, DR);
        for (int i=0; i<d; i++)
        {
            for (int j=0; j<DL; j++)
            {
                for (int k=0; k<DR; k++)
                {
                    fold (i*DL+j,k) = tsrA.getTensor(j,k,i);
                }
            }
        }

        /**
         * Singular Value Decomposition 
         */
        DB = d*DL > DR ? DR : d*DL;
        JacobiSVD<MatrixXcd> svd (fold, ComputeThinU | ComputeThinV);

        //  demision consistency test
        MatrixXcd U = svd.matrixU ();
        MatrixXcd V = svd.matrixV ().transpose ();    //USV'

        //  Calculate S=S*V
        S.resize (DB, DB);  
        S = MatrixXcd::Zero (DB,DB);
        for (int i=0; i<DB; i++) S (i,i) = svd.singularValues ()(i,0);
        S = S*V;

        /**
         * write U to tensor projector tsrA
         */
        tsrA.reAlloc (DL, DB);
        for (int k=0; k<d; k++)
        {
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DB; j++)
                {
                    tsrA.setTensor (i, j, k, U (k*DL+i,j));
                }
            }
        }

        return;
    }
    else
    {
        //  prepare from right to left

        /**
         * construct matrix fold with dim (DL, d*DR)
         */
        MatrixXcd fold = MatrixXcd::Zero (DL, d*DR);
        for (int j=0; j<DL; j++)
        {
            for (int i=0; i<d; i++)
            {
                for (int k=0; k<DR; k++)
                {
                    fold (j,i*DR+k) = tsrA.getTensor (j,k,i);
                }
            }
        }

        /**
         * Singular Value Decomposition 
         */
        DB = DL > d*DR ? d*DR : DL;
        JacobiSVD<MatrixXcd> svd (fold, ComputeThinU | ComputeThinV);

        //  demision consistency test
        MatrixXcd U = svd.matrixU ();
        MatrixXcd V = svd.matrixV ().transpose ();    //USV'
      //  cout << "\nDB = " << DB << " | V = \n" << V << endl;

        //  Calculate S=U*S
        S.resize (DB, DB);  
        S = MatrixXcd::Zero (DB,DB);
        for (int i=0; i<DB; i++) S (i,i) = svd.singularValues ()(i,0);
        S = U*S;

        /**
         * write V to tensor projector tsrA 
         */
        tsrA.reAlloc (DB, DR);
        for (int k=0; k<d; k++)
        {
            for (int i=0; i<DB; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    tsrA.setTensor (i, j, k, V (i, k*DR+j));
                }
            }
        }

        return;
    }

    //  end prepare_onesite
    return;
}

void VMPS::storage_update (
        const bool& dir,    //  0=lr, 1=rl
        const int& sj,
        const TensorProj& tsrA,
        const bool& gs,
        const TensorProj& tsrB
        )
{
    const int DLA = tsrA.getDL ();
    const int DRA = tsrA.getDR ();
    const int dA  = tsrA.getd  ();

    const int DLB = tsrB.getDL ();
    const int DRB = tsrB.getDR ();
    const int dB  = tsrB.getd  ();

    /**
     * store matrix A (DLA*DRA, dA) corresponding to tsrA
     * store matrix B (DLB*DRB, dB) corresponding to tsrB
     */
    MatrixXcd A  = MatrixXcd::Zero (DLA*DRA, dA);
    for (int i=0; i< DLA; i++)
    {
        for (int j=0; j< DRA; j++)
        {
            for (int k=0; k<dA; k++)
            {
                A (i*DRB+j,k) = tsrA.getTensor (i,j,k);
            }
        }
    }

    MatrixXcd B  = MatrixXcd::Zero (DLB*DRB, dB);
    if (gs)
    {
        for (int i=0; i< DLB; i++)
        {
            for (int j=0; j< DRB; j++)
            {
                for (int k=0; k<dB; k++)
                {
                    B (i*DRB+j,k) = tsrB.getTensor (i,j,k);
                }
            }
        }
    }

    //  update 

    if (dir == false)   //lr
    {
        /*****************
         * Update Hstorage
         * ***************/
        for (int m=0; m<M; m++)
        {
            /**
             * contract A and h and A+
             */
            MatrixXcd U1 = A * hset[m][sj] * A.adjoint ();

            /**
             * reshape U1 to have dim (DLA*DLA, DRA*DRA)
             */
            MatrixXcd S1= MatrixXcd::Zero (DLA*DLA, DRA*DRA);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLA; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRA; l++)
                        {
                            S1(i*DLA+j, k*DRA+l) = U1(i*DRA+k, j*DRA+l);
                        }
                    }
                }
            }

            /**
             * contract S1 with Hstorage[m][sj]
             */
            Hstorage[m][sj+1] = Hstorage[m][sj]*S1;

        }// end m

        /*****************
         * Update Nstorage
         * ***************/

        if (pbc || !USE_PREPARE)
        {
            /**
             * contract A and h and A+
             */
            MatrixXcd U2= A * A.adjoint ();

            /**
             * reshape U2 to have dim (DLA*DLA, DRA*DRA)
             */
            MatrixXcd S2= MatrixXcd::Zero (DLA*DLA, DRA*DRA);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLA; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRA; l++)
                        {
                            S2 (i*DLA+j, k*DRA+l) = U2 (i*DRA+k, j*DRA+l);
                        }
                    }
                }
            }

            /**
             * contract S2 with Nstorage[sj]
             */
            Nstorage[sj+1] = Nstorage[sj]*S2;
        }

        /*****************
         * Update Cstorage
         * ***************/

        if (gs)
        {
            /**
             * contract A and B -- <B|A>
             */
            MatrixXcd U3 = A*B.adjoint ();

            //  (3) reshape U3 to have dim (DLA*DLB, DRA*DRB)
            MatrixXcd S3 = MatrixXcd::Zero (DLA*DLB, DRA*DRB);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLB; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRB; l++)
                        {
                            S3 (i*DLB+j, k*DRB+l) = U3 (i*DRA+k, j*DRB+l);
                        }
                    }
                }
            }

            //  contract S3 with Cstorage[sj]
            Cstorage[sj+1] = Cstorage[sj]*S3;
        }

        return;
    }
    else    //rl
    {
        /*****************
         * Update Hstorage
         * ***************/
        for (int m=0; m<M; m++)
        {
            /**
             * contract A and h and A+
             */
            MatrixXcd U1 = A * hset[m][sj] * A.adjoint ();

            /**
             * reshape U1 to have dim (DLA*DLA, DRA*DRA)
             */
            MatrixXcd S1 = MatrixXcd::Zero (DLA*DLA, DRA*DRA);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLA; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRA; l++)
                        {
                            S1 (i*DLA+j, k*DRA+l) = U1 (i*DRA+k, j*DRA+l);
                        }
                    }
                }
            }

            //  (4) contract S1 with Hstorage[m][sj+1]
            Hstorage[m][sj] = S1 * Hstorage[m][sj+1];
        }

        /*****************
         * Update Nstorage
         * ***************/

        if (pbc || !USE_PREPARE)
        {
            /**
             * contract A and h and A+
             */
            MatrixXcd U2 = A * A.adjoint ();

            /**
             * reshape U2 to have dim (DLA*DLA, DRA*DRA)
             */
            MatrixXcd S2 = MatrixXcd::Zero (DLA*DLA, DRA*DRA);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLA; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRA; l++)
                        {
                            S2 (i*DLA+j, k*DRA+l) = U2 (i*DRA+k, j*DRA+l);
                        }
                    }
                }
            }

            /**
             * contract S2 with Nstorage[sj+1]
             */
            Nstorage[sj] = S2 * Nstorage[sj+1];
        }

        /*****************
         * Update Cstorage
         * ***************/
        if (gs)
        {
            /**
             * contract A and B -- <B|A>
             */
            MatrixXcd U3 = A*B.adjoint ();

            //  (3) reshape U3 to have dim (DLA*DLB, DRA*DRB)
            MatrixXcd S3 = MatrixXcd::Zero (DLA*DLB, DRA*DRB);
            for (int i=0; i<DLA; i++)
            {
                for (int j=0; j<DLB; j++)
                {
                    for (int k=0; k<DRA; k++)
                    {
                        for (int l=0; l<DRB; l++)
                        {
                            S3 (i*DLB+j, k*DRB+l) = U3 (i*DRA+k, j*DRB+l);
                        }
                    }
                }
            }

            //  contract S3 with Cstorage[sj+1]
            Cstorage[sj] = S3 * Cstorage[sj+1];
        }

        return;
    }

    //end
    return;
}

double VMPS::Mxyz (const char& ch, const bool& pbc) const
{
    //  choose pauli matrix
    MatrixXcd pauli (2,2);
    if (ch == 'x') pauli = sx;
    else if (ch == 'y') pauli = sy;
    else if (ch == 'z') pauli = sz;
    else pauli = id;

    //  result
    double mag = 0.0;

    //  contraction from left to right
    for (int n1=0; n1<N; n1++)
    {
        //  calc <psi|sigma(i)|psi>
        int DL0 = (*mpsGS).getTsrProj (0).getDL ();
        MatrixXcd E = MatrixXcd::Identity (DL0*DL0, DL0*DL0);

        for (int n2=0; n2<N; n2++)
        {
            const TensorProj& tsr = (*mpsGS).getTsrProj (n2);
            const int DL = tsr.getDL ();
            const int DR = tsr.getDR ();
            const int d  = tsr.getd  ();

            MatrixXcd U = MatrixXcd::Zero (DL*DR,d);
            //  contract U with pauli and U+
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        U (i*DR+j, k) = tsr.getTensor (i,j,k); 
                    }
                }       
            }

            //  contract U with pauli matrix and U+
            MatrixXcd V (DL*DR, DL*DR);
            if (n1 == n2) 
            {
                /*
                cout << "\nU.cols () = " << U.cols () << 
                        "  Pauli.rows () = " << pauli.rows () << 
                        "  Pauli.cols () = " << pauli.cols () << endl;

                assert (U.cols () == pauli.rows ());
                assert (U.cols () == pauli.cols ()); 
                */

                V = U * pauli * U.adjoint ();
            }
            else V = U*U.adjoint ();

            //  reshape V to P
            MatrixXcd P = MatrixXcd::Zero (DL*DL,DR*DR);
            for (int i1=0; i1<DL; i1++)
            {
                for (int j1=0; j1<DR; j1++)
                {
                    for (int i2=0; i2<DL; i2++)
                    {
                        for (int j2=0; j2<DR; j2++)
                        {
                            P (i1*DL+i2, j1*DR+j2) = V (i1*DR+j1, i2*DR+j2);
                        }
                    }
                }
            }

            //  contract
            /*
            cout << "\nE.cols () = " << E.cols () << " P.rows () = " << P.rows () << endl;
            assert (E.cols () == P.rows ());
            */

            E = E*P;
        }

        mag += E.trace ().real ();
    }

    return mag/N;
}

double VMPS::Cxyz (const char& ch, const int& r, const bool& pbc) const
{
    //  choose pauli matrix
    MatrixXcd pauli (2,2);
    if (ch == 'x') pauli = sx;
    else if (ch == 'y') pauli = sy;
    else if (ch == 'z') pauli = sz;
    else pauli = id;

    //  result
    double corr = 0.0;

    //  contraction from left to right
    for (int n1=0; n1<N; n1++)
    {
        //  calc <psi|sigma(i)|psi>
        int DL0 = (*mpsGS).getTsrProj (0).getDL ();
        MatrixXcd E = MatrixXcd::Identity (DL0*DL0, DL0*DL0);

        //  n2^{prime}
        int n2p = 0;
        if (!pbc && n1+r >= N) continue;
        else n2p = (n1+r) % N;

        for (int n2=0; n2<N; n2++)
        {
            const TensorProj& tsr = (*mpsGS).getTsrProj (n2);
            const int DL = tsr.getDL ();
            const int DR = tsr.getDR ();
            const int d  = tsr.getd  ();

            MatrixXcd U = MatrixXcd::Zero (DL*DR,d);
            //  contract U with pauli and U+
            for (int i=0; i<DL; i++)
            {
                for (int j=0; j<DR; j++)
                {
                    for (int k=0; k<d; k++)
                    {
                        U (i*DR+j, k) = tsr.getTensor (i,j,k); 
                    }
                }       
            }

            //  contract U with pauli matrix and U+
            MatrixXcd V (DL*DR, DL*DR);
            if (n2==n1 || n2==n2p) V = U*pauli*U.adjoint ();
            else V = U*U.adjoint ();

            //  reshape V to P
            MatrixXcd P = MatrixXcd::Zero (DL*DL,DR*DR);
            for (int i1=0; i1<DL; i1++)
            {
                for (int j1=0; j1<DR; j1++)
                {
                    for (int i2=0; i2<DL; i2++)
                    {
                        for (int j2=0; j2<DR; j2++)
                        {
                            P (i1*DL+i2, j1*DR+j2) = V (i1*DR+j1, i2*DR+j2);
                        }
                    }
                }
            }

            //  contract
            E = E*P;
        }

        corr += E.trace ().real ();
    }

    return corr/N;
}

int VMPS::Output (const RI& ri)
{
    time_t start, end;
    time (&start);

    //  ofstream
    //  file name must be "output_res.qs"
    ofstream outfile;
    outfile.open ("output_res.qs");
    if (!outfile)
    {
        cerr << "error: unable to open output file: "
             << outfile << endl;
        return -1;
    }

    //  currently unavailable
    outfile << setiosflags (ios::fixed) << setprecision (10);

    outfile << "******************************************************";
    outfile << "\nYOU ARE USING PROGRAM: qsAnnni -- written by Chen Liao" << endl;
    outfile << "\nReview of Input Variables:" << endl
            << "\nNS =   " << ri.NS
            << "\nVD    =   " << ri.VD 
            << "\nCL    =   " << ri.CL
            << "\nJ1X   =   " << ri.J1X
            << "\nJ1Y   =   " << ri.J1Y
            << "\nJ1Z   =   " << ri.J1Z
            << "\nJ2X   =   " << ri.J2X
            << "\nJ2Y   =   " << ri.J2Y
            << "\nJ2Z   =   " << ri.J2Z
            << "\nHX    =   " << ri.HX
            << "\nHY    =   " << ri.HY
            << "\nHZ    =   " << ri.HZ
            << endl;

    //  output energy
    outfile << "\ngroud state energy is: " << E0 << endl;

    //  output energy gap
    if (ri.CALC_FES)
    {
        double gap = E1-E0;
        outfile << "\nfirst excited state energy is: " << E1
                << "\nenergy gap between gs and fes is: " << gap << endl;
    }

    //  output magnetization
    double Mx, My, Mz;
    Mx = My = Mz = 0.0;
    if (ri.OUT_MAG) 
    {
        Mx = Mxyz ('x', ri.USE_PBC);
        My = Mxyz ('y', ri.USE_PBC);
        Mz = Mxyz ('z', ri.USE_PBC);

        outfile << "\nmagnetization along X axis is: " << Mx
                << "\nmagnetization along Y axis is: " << My
                << "\nmagnetization along Z axis is: " << Mz
                << endl;
    }

    //  output correlation
    double Cx, Cy, Cz;
    Cx = Cy = Cz = 0.0;
    if (ri.OUT_CORR)  
    {
        Cx = Cxyz ('x', ri.CL, ri.USE_PBC);
        Cy = Cxyz ('y', ri.CL, ri.USE_PBC);
        Cz = Cxyz ('z', ri.CL, ri.USE_PBC);

        outfile << "\ncorrelation (length " << ri.CL 
                << ") along X axis is: " << Cx
                << "\ncorrelation (length " << ri.CL 
                << ") along Y axis is: " << Cy
                << "\ncorrelation (length " << ri.CL
                << ") along Z axis is: " << Cz
                << endl;
    }

    //  output all data in one line
    outfile << "\nFormat:\n" << E0 << " " << E1 
            << " " << Mx << " " << My << " " << Mz
            << " " << Cx << " " << Cy << " " << Cz 
            << endl;

    //  output cfgs
    
    CFGSTR amp;
    getConfig (*mpsGS, amp);
    CFGSTR::const_reverse_iterator it = amp.rbegin ();

    double prob = 0.0;
    int lc = 0;
    while (it != amp.rend ())
    {
        outfile << "\nbasis " << setw(2) << ++lc << ", amp = " << setw(24) << it->second
                << ", prob = " << norm (it->second) << " | " << it->first;
        prob += norm (it->second); 
        it++;
    }

    assert (fabs (1.0-prob) < TINY);
    cout << endl;

    /**
     * the section of codes are for test : analytical calculation of Mx
     */
    
    //  test begin
    /*
    double testMx = 0.0;
    double testMz = 0.0;

    it = amp.rbegin ();
    while (it != amp.rend ())
    {
        string cfgstr (it->first);
        complex<double> coef = it->second;

        int signMz = 0;
        for (int j=0; j<N; j++)
        {
            //  Mx
            string cfgstr_conj (cfgstr);
            cfgstr_conj[j] = (cfgstr[j] == 'u' ? 'd':'u');

            CFGSTR::const_iterator it_find = amp.find (cfgstr_conj);
            assert (it_find != amp.end ());

            testMx += (conj (coef) * it_find->second).real ();

            //  Mz
            signMz += (cfgstr[j] == 'u' ? 1:-1);
        }
        
        testMz += (conj (coef) * coef).real () * signMz;

        it++;
    }

    outfile << "\n\ntest: Mx = " << testMx / N 
            << "\ntest: Mz = " << testMz / N << endl;
    */
    //  test end

    //  time output
    time (&end);
    timeOUT = difftime (end,start);
    
    outfile << endl << endl;
    outfile << "\nTime Consumption:" << endl
            << "\nVMPS::Opt: " << timeOPT  << "s"
            << "\nDiag::OptOne: "   << timeOPT_one  << "s"
            << "\nDiag::Output: "  << timeOUT << "s"
            << endl << endl;

    outfile.close ();

    return 0;
}

void VMPS::getConfig (const MPS& mpsRef, CFGSTR& cfgmap) const
{
    /**
     * output configurations
     * array script: 0 --> up, 1 --> down
     */

    //  allocation
    int cfg_tot = pow (2, N);
    char* config = new char[N];
    for (int i=0; i<N; i++) config[i] = '0'; // all up, FM

    //  calculation
    for (int i=0; i<cfg_tot; i++)
    {
        int DL0 = mpsRef.getTsrProj (0).getDL ();
        MatrixXcd E = MatrixXcd::Identity (DL0, DL0);

        for (int j=0; j<N; j++)
        {
            const TensorProj& tsr = (*mpsGS).getTsrProj (j);
            const int DL = tsr.getDL ();
            const int DR = tsr.getDR ();
            const int d  = tsr.getd  ();

            MatrixXcd U = MatrixXcd::Zero (DL,DR);
            int ind_d = (config[j] == '0'? 0:1);

            for (int k=0; k<DL; k++)
            {
                for (int l=0; l<DR; l++)
                {
                    U (k, l) = tsr.getTensor (k,l,ind_d); 
                }       
            }

            E = E*U;
        }

        complex<double> amp = E.trace ();

        //  '0' --> up, '1' --> down
        string cfgAstr ("");
        for (int j=0; j<N; j++) cfgAstr += (config[j] == '0' ? 'u':'d');

        //  '0' --> up, '1' --> down
        cfgmap.insert (CFGSTR::value_type (cfgAstr,amp));

        //  find new configurations
        int count = N-1;
        bool iflag = 1;

        while (iflag && count >= 0) 
        {
            if (config[count] == '0') {config[count] = '1'; iflag = 0;}
            else {config[count] = '0'; count --; iflag = 1;}
        }
    }

    delete [] config;

    return;
}
