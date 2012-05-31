#include "diag.h"

Diag::Diag ()
{
    sx = sy = sz = id = MatrixXcd::Zero (2,2);
    sx (1,0) = sx (0,1) = 1.0;
    sy (0,1) = complex<double>(0.0,-1.0); sy (1,0) = complex<double>(0.0,1.0);
    sz (0,0) = 1.0; sz (1,1) = -1.0;
    id (0,0) = id (1,1) = 1.0;
    timeMH = timeDH = timeOUT = 0.0;
}

Diag::~Diag () {} 

int Diag::Output (const RI& ri)
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

    //  output eigenvalues/eigenvectors first
 //   outfile << "\nThe eigenvalues of H are:\n"  << eigenval << endl;
 //   outfile << "\nThe eigenvectors of H are:\n" << eigenvec << endl;
    
    outfile << setiosflags (ios::fixed) << setprecision (10);

    //  configuration output

    //  generate configuration notations
    vector<string> cfgTmp1, cfgTmp2;
    cfgTmp1.push_back ("u");
    cfgTmp1.push_back ("d");

    int currD = 2; 

    for (int i=1; i<ri.NS; i++)
    {
        for (int j=0; j<currD; j++) 
        {
            string str = cfgTmp1[j];
            cfgTmp2.push_back (str + 'u');
            cfgTmp2.push_back (str + 'd');
        }

        cfgTmp1.clear ();
        cfgTmp1 = cfgTmp2;
        cfgTmp2.clear ();

        currD *= 2;
    }

    outfile << "\nConfigurations are:" << endl;
    for (int i=0; i<eigenval.rows (); i++)
    {
        outfile << "\nLevel = " << i << " | Energy = " << eigenval (i);
        outfile << "\n------------------------";

   //     double prob_tot = 0.0;
        for (int j=0; j<eigenval.rows (); j++)
        {
            double prob = norm (eigenvec (j,i));
            if (prob < TINY_Q)
            {
                eigenvec (j,i).real () = eigenvec (j,i).imag () = 0.0;
                prob = 0.0;
            }

     //       prob_tot += prob;
            outfile << "\nbasis " << setw (2) << j+1 << ", amp = " << setw(22) 
                    <<  eigenvec (j,i) << ", prob = " << prob << " | " << cfgTmp1[j];
        }

   //     cout << "\nLevel = " << i << " prob in total = " << prob_tot << endl;
     //   assert (fabs(1.0-prob_tot) < TINY_Q);
        outfile << endl;
    }


    //  output energy gap
    double gap = eigenval (1)-eigenval (0);
    outfile << "\nenergy gap between gs and fes is: " << gap << endl;

    //  output magnetization
    if (ri.OUT_MAG)
    {
        outfile << "\nmagnetization along X axis is: " << Mxyz('x', ri.NS, 0) << setw (14) << Mxyz('x', ri.NS, 1)
                << "\nmagnetization along Y axis is: " << Mxyz('y', ri.NS, 0) << setw(14) << Mxyz('y', ri.NS, 1)
                << "\nmagnetization along Z axis is: " << Mxyz('z', ri.NS, 0) << setw(14) << Mxyz('z', ri.NS, 1)
                << endl;
    }

    //  output correlation
    if (ri.OUT_CORR)
    {
         outfile << "\ncorrelation (length " << ri.CL 
                 << ") along X axis is: " << Cxyz('x', ri.CL, ri.NS, 0, ri.USE_PBC) << setw(14) << Cxyz('x', ri.CL, ri.NS, 1, ri.USE_PBC)
                 << "\ncorrelation (length " << ri.CL 
                 << ") along Y axis is: " << Cxyz('y', ri.CL, ri.NS, 0, ri.USE_PBC) << setw(14) << Cxyz('y', ri.CL, ri.NS, 1, ri.USE_PBC)
                 << "\ncorrelation (length " << ri.CL 
                 << ") along Z axis is: " << Cxyz('z', ri.CL, ri.NS, 0, ri.USE_PBC) << setw(14) << Cxyz('z', ri.CL, ri.NS, 1, ri.USE_PBC)
                 << endl;
    }

    //  time output
    time (&end);
    timeOUT = difftime (end,start);
    
    double timeTOT = timeMH + timeDH + timeOUT;
    outfile << "\nDiag::MatrixH: " << timeMH  << "s"
            << "\nDiag::DiagH: "   << timeDH  << "s"
            << "\nDiag::Output: "  << timeOUT << "s"
            << "\nTime in total for class Diag: " << timeTOT << "s"
            << endl;

    outfile.close ();
    return 0;
}

double Diag::Mxyz (const char& xyz, const int& N, const int& level)
{
    //  result
    double res = 0.0;

    //  x, y, z
    MatrixXcd sxyz (2,2);
    if (xyz == 'x') sxyz = sx;
    else if (xyz == 'y') sxyz = sy;
    else if (xyz == 'z') sxyz = sz;
    else sxyz = id;
    
    //  eigenvec 
    VectorXcd vec = eigenvec.col (level);
 //   VectorXcd vec = 1.0/sqrt(2.0)*(eigenvec.col(0)-eigenvec.col(1));

    //  magnetization
    for (int i=0; i<N; i++)
    {
        //  expected value of each s(x/y/z)

        MatrixXcd Xtmp2;
        if (i == 0) Xtmp2 = sxyz;
        else Xtmp2 = id;

        int currD = 2;

        //  1*1...s(i)*1*1...
        for (int j=1; j<N; j++)
        {
            currD *= 2;
            MatrixXcd Xtmp1 = MatrixXcd::Zero (currD, currD);

            if (j == i) cross (Xtmp2, sxyz, Xtmp1);
            else cross (Xtmp2, id, Xtmp1);

            Xtmp2 = Xtmp1;
        }
        
        // calculate expected value
        complex<double> eval = vec.adjoint () * Xtmp2 * vec;
        res += eval.real ();
    }

    return res/N;
}

double Diag::Cxyz (const char& xyz, const int& dis, const int& N, const int& level, const bool& pbc)
{
    //  result
    double res = 0.0;

    //  x, y, z
    MatrixXcd sxyz (2,2);
    if (xyz == 'x') sxyz = sx;
    else if (xyz == 'y') sxyz = sy;
    else if (xyz == 'z') sxyz = sz;
    else sxyz = id;

    //  eigenvec 
    VectorXcd vec = eigenvec.col (level);

    //  correlation
    for (int i=0; i<N; i++)
    {
        //  break for OBC
        if (i+dis > N-1 && !pbc) break;

        //  interaction site index
        const int idis = (i+dis) % N;

        MatrixXcd Xtmp2;
        if (i == 0 || idis == 0) Xtmp2 = sxyz;
        else Xtmp2 = id;

        int currD = 2;

        //  1*1...s(i)*..1..s(idis)...*1...
        for (int j=1; j<N; j++)
        {
            currD *= 2;
            MatrixXcd Xtmp1 = MatrixXcd::Zero (currD, currD);

            if (j == i || j == idis) cross (Xtmp2, sxyz, Xtmp1);
            else cross (Xtmp2, id, Xtmp1);

            Xtmp2 = Xtmp1;
        }
        
        // calculate expected value
        complex<double> eval = vec.adjoint () * Xtmp2 * vec;
        res += eval.real ();
    }

    return res/N;
}


void Diag::DiagH ()
{
    //  time
    time_t start, end;
    time (&start);

    //  eigensolver
    SelfAdjointEigenSolver<MatrixXcd> eigensolver(H);
    eigenval = eigensolver.eigenvalues ();
    eigenvec = eigensolver.eigenvectors ();

    //  print
    cout << "\nenergy of ground state : " << setiosflags (ios::fixed) 
         << setprecision (10) << eigenval (0,0);
    cout << "\nenergy of first excited state : " << setiosflags (ios::fixed)
         << setprecision (10) << eigenval (1,0) << "\n\n" << endl;
//    cout << "\nThe eigenvalues of H are:\n"  << eigenval << endl;
//  cout << "\nThe eigenvectors of H are:\n" << eigenvec << endl;

    //  time
    time (&end);
    timeDH = difftime (end,start);
    
    return;
}

void Diag::MatrixH (const RI& ri)
{
    //  time
    time_t start, end;
    time (&start);

    //  initialization of H
    const int dimH = pow (2, ri.NS);
    H = MatrixXcd::Zero (dimH, dimH);

    //  calculate elements of H
    if (fabs (ri.J1X) > TINY_Q) H_J ('x', 1, ri.NS, ri.J1X, ri.USE_PBC); //J1X
    if (fabs (ri.J1Y) > TINY_Q) H_J ('y', 1, ri.NS, ri.J1Y, ri.USE_PBC); //J1Y
    if (fabs (ri.J1Z) > TINY_Q) H_J ('z', 1, ri.NS, ri.J1Z, ri.USE_PBC); //J1Z
    if (fabs (ri.J2X) > TINY_Q) H_J ('x', 2, ri.NS, ri.J2X, ri.USE_PBC); //J2X
    if (fabs (ri.J2Y) > TINY_Q) H_J ('y', 2, ri.NS, ri.J2Y, ri.USE_PBC); //J2X
    if (fabs (ri.J2Z) > TINY_Q) H_J ('z', 2, ri.NS, ri.J2Z, ri.USE_PBC); //J2X
    if (fabs (ri.HX)  > TINY_Q) H_M ('x', ri.NS, ri.HX); //Hx
    if (fabs (ri.HY)  > TINY_Q) H_M ('y', ri.NS, ri.HY); //Hy
    if (fabs (ri.HZ)  > TINY_Q) H_M ('z', ri.NS, ri.HZ); //Hz

    //  print H
 //   cout << "\nHamiltonian Matrix :\n" << H << endl;

    //  time
    time (&end);
    timeMH = difftime (end,start);
    
    return;
}

void Diag::H_J (
        const char& xyz, 
        const int&  dis, 
        const int&  N, 
        const double& Jdis, 
        const bool& pbc
        )
{
    //  J(X,Y,Z)
    MatrixXcd sxyz (2,2);
    if (xyz == 'x') sxyz = sx;
    else if (xyz == 'y') sxyz = sy;
    else if (xyz == 'z') sxyz = sz;
    else sxyz = id;

    //  H(j1) = -J \sum_{i} s(i)*s(i+1)
    for (int i=0; i<N; i++)
    {
        //  break for OBC
        if (i+dis > N-1 && !pbc) break;

        //  interaction site index
        const int idis = (i+dis) % N;

        MatrixXcd Xtmp2;
        if (i == 0 || idis == 0) Xtmp2 = sxyz;
        else Xtmp2 = id;

        int currD = 2;

        //  1*1...s(i)*..1..s(idis)...*1...
        for (int j=1; j<N; j++)
        {
            currD *= 2;
            MatrixXcd Xtmp1 = MatrixXcd::Zero (currD, currD);

            if (j == i || j == idis) cross (Xtmp2, sxyz, Xtmp1);
            else cross (Xtmp2, id, Xtmp1);

            Xtmp2 = Xtmp1;
        }
        
        //  write to H
        H -= Jdis*Xtmp2;
    }

    return;
}

void Diag::H_M (const char& xyz, const int& N, const double& hxyz)
{
    //  H (X,Y,Z)
    MatrixXcd sxyz (2,2);
    if (xyz == 'x') sxyz = sx;
    else if (xyz == 'y') sxyz = sy;
    else if (xyz == 'z') sxyz = sz;
    else sxyz = id;

    //  H(m) = -h * \sum_{i} s(i)
    for (int i=0; i<N; i++)
    {
        MatrixXcd Xtmp2;
        if (i == 0) Xtmp2 = sxyz;
        else Xtmp2 = id;

        int currD = 2;

        //  1*1...s(i)*1*1...
        for (int j=1; j<N; j++)
        {
            currD *= 2;
            MatrixXcd Xtmp1 = MatrixXcd::Zero (currD, currD);

            if (j == i) cross (Xtmp2, sxyz, Xtmp1);
            else cross (Xtmp2, id, Xtmp1);

            Xtmp2 = Xtmp1;
        }
        
        //  write to H
        H -= hxyz * Xtmp2;
    }

    return;
}

void Diag::cross (const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C) const
{
    //  C must be initialized with ZERO
    //  C(m*p,n*q) = A(m,n) product B(p,q)

    //  Brow,Bcol
    const int& Brow = B.rows ();
    const int& Bcol = B.cols ();
    const int& Bmin = Brow>Bcol ? Bcol:Brow;

    //  whether B is equal to Identity
    bool BeqI = B.isIdentity ();
    if (BeqI)
    {
        for (int i=0; i<A.rows (); i++)
        {
            for (int j=0; j<A.cols (); j++)
            {
                for (int k=0; k<Bmin; k++)
                {
                    C(i*Brow+k,j*Bcol+k) = A (i,j);
                }
            }
        }
    }
    else
    {
        //  block operation
        for (int i=0; i<A.rows (); i++)
        {
            for (int j=0; j<A.cols (); j++)
            {
                for (int k=0; k<Brow; k++)
                {
                    for (int l=0; l<Bcol; l++)
                    {
                        C(i*Brow+k, j*Bcol+l) = A(i,j)*B(k,l);
                    }
                }
            }
        }
    }

    return;
}
