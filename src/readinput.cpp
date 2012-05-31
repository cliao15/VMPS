#include "readinput.h"

int RI::readInputVar ()
{
    //  file name must be "input.qs"
    ifstream infile;
    infile.open ("input.qs");
    if (!infile)
    {
        cerr << "error: unable to open input file: "
             << infile << endl;
        return -1;
    }

    //  streamsize
    const int strms = 256;

    //  lines starting with '#' are considered as a comment and ignored
    //  however, lines containing '#' not as the first char would not be igonred
    char ch;
    while (infile >> ch, ! infile.eof ())
    {
        //  test
  //      cout << "\nThe first letter is " << ch;
        
        //  stream in good state?
        if (infile.bad ()) throw runtime_error ("IO stream corrupted");
        if (infile.fail ()) throw runtime_error ("bad data");
        if (ch == '#')
        {
            //  read this line and ignore it
            char comment[256];
            infile.getline (comment,256);

            //  test
      //      cout << "\nComment ignored is " << comment;
        }
        else
        {
            //  read input variable
            infile.unget ();    //set the last read char available again for the next read
            char data[256];
            infile.getline (data,256,' ');

            //  clear spaces
            infile >> ch;
            while (ch == ' ') infile >> ch;
            infile.unget ();
            
            //  read val
            if (string (data) == "NS")       infile >> NS;
            if (string (data) == "VD")          infile >> VD;
            if (string (data) == "MTH")      infile >> MTH;
            if (string (data) == "INIT_CFG")         infile >> INIT_CFG;
            if (string (data) == "CL")  infile >> CL;
            if (string (data) == "J1X")         infile >> J1X;
            if (string (data) == "J1Y")         infile >> J1Y;
            if (string (data) == "J1Z")         infile >> J1Z;
            if (string (data) == "J2X")         infile >> J2X;
            if (string (data) == "J2Y")         infile >> J2Y;
            if (string (data) == "J2Z")         infile >> J2Z;
            if (string (data) == "HX")          infile >> HX;
            if (string (data) == "HY")          infile >> HY;
            if (string (data) == "HZ")          infile >> HZ;
            if (string (data) == "TOL")         infile >> TOL;
            if (string (data) == "CALC_FES")         infile >> CALC_FES;
            if (string (data) == "OUT_MAG")        infile >> OUT_MAG;
            if (string (data) == "OUT_CORR")         infile >> OUT_CORR;
            if (string (data) == "OUT_CFG")         infile >> OUT_CFG;
            if (string (data) == "USE_PBC")         infile >> USE_PBC;
        }
    }

    //  close
    infile.close ();

    printVal ();
    return 0;
}

int RI::printVal ()
{
    //  file name must be "output_val.qs"
    ofstream outfile;
    outfile.open ("output_val.qs");
    if (!outfile)
    {
        cerr << "error: unable to open output file: "
             << outfile << endl;
        return -1;
    }

    //  printing
    outfile << "Review of Input Variables:" << endl
            << "\nNS =   " << NS
            << "\nVD    =   " << VD 
            << "\nMTH    =   " << MTH
            << "\nINIT_CFG   =   " << INIT_CFG
            << "\nCL    =   " << CL
            << "\nJ1X   =   " << J1X
            << "\nJ1Y   =   " << J1Y
            << "\nJ1Z   =   " << J1Z
            << "\nJ2X   =   " << J2X
            << "\nJ2Y   =   " << J2Y
            << "\nJ2Z   =   " << J2Z
            << "\nHX    =   " << HX
            << "\nHY    =   " << HY
            << "\nHZ    =   " << HZ
            << "\nTOL   =   " << TOL
            << "\nCALC_FES   =   " << CALC_FES
            << "\nUSE_PBC   =   " << USE_PBC
            << "\nOUT_MAG  =   " << OUT_MAG
            << "\nOUT_CORR   =   " << OUT_CORR
            << "\nOUT_CFG   =   " << OUT_CFG
            << endl;
    outfile.close ();
    return 0;
}
