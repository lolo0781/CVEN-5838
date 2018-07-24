//  ================================================================================
//  ||                                                                            ||
//  ||              fd                                                            ||
//  ||              -------------------------------------------                   ||
//  ||              F I N I T E   D I F F E R E N C E                             ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              -------------------------------------------                   ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     Los Alamos National Laboratory                         ||
//  ||                                                                            ||
//  ||                For: CU-Boulder CVEN 5838-200 and -200B                     ||
//  ||                     (Not for distribution outside of class)                ||
//  ||                                                                            ||
//  ||          Copyright: 2017-2018 Scott Runnels                                ||
//  ||                                                                            ||
//  ================================================================================

#include "fd.h"

//  ==
//  ||
//  ||      C L A S S:   L A P L A C I A N O N G R I D 
//  ||
//  ||      Does the following:
//  ||
//  ||      (1) Forms a linear system Ax=b representing a 2D
//  ||          finite difference approximation to Laplaces
//  ||          equation on a square.
//  ||      (2) Performs Gauss-Seidel iteration on that linear
//  ||          system up to a specified number of iterations
//  ||          or a specified tolerance.
//  ||          
//  ||       Here is what the grid looks like:
//  ||          
//  ||          
//  ||                                                     ny*nx
//  ||       *-------*-------*-------*-------*  ...   -------*
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||       .       .       .       .       .               .
//  ||       .       .       .       .       .               .
//  ||       .       .       .       .       .               .
//  || 2nx+1  
//  ||       *-------*-------*-------*-------*  ...   -------*
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||  nx+1 |       |       |       |       |          2nx  |
//  ||       *-------*-------*-------*-------*  ...   -------*
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||       |       |       |       |       |               |
//  ||       *-------*-------*-------*-------*  ...   -------*
//  ||      1       2       3       4                       nx  
//  ||    
//  == 

class LaplacianOnGrid
{

public:

  int nx    , ny   , nrows;
  _D_ length, dx;
  VDD A ;    VD  phi ;  VD  b ;
  string Name;

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  LaplacianOnGrid(int ncell_x, int ncell_y, string _Name )
  {
    nx     = ncell_x + 1;
    ny     = ncell_y + 1;
    nrows  = nx*ny;
    length = 1.;
    dx     = length / (nx-1);
    Name   = _Name;

    A.resize(nrows+1); rLOOP A[r].resize(nrows+1);
    b.resize(nrows+1);
    phi.resize(nrows+1);
  }

  //  ==
  //  ||
  //  ||  Form Linear System Ax = b
  //  ||
  //  ==

  void FormLS()
  {
    
    // ---------------------------------------------------
    // Initialize linear system
    // ---------------------------------------------------
    
    rLOOP cLOOP A[r][c] = 0.;
    rLOOP b[r] = 0.;
    rLOOP phi[r] = 0.;
    
    // ---------------------------------------------------
    // Form matrix entries for the interior grid points
    // ---------------------------------------------------

    _D_ dx2 = dx*dx;
    
    iLOOPi 
      jLOOPi
      {
	int p = pid(i,j);
	A[ p ][ p    ] = -4./dx2;
	A[ p ][ p+1  ] =  1./dx2;
	A[ p ][ p-1  ] =  1./dx2;
	A[ p ][ p+nx ] =  1./dx2;
	A[ p ][ p-nx ] =  1./dx2;
      }
      
    // ---------------------------------------------------
    // Apply boundary conditions
    // ---------------------------------------------------
    
    iLOOP
      {
	int p,j;
	p = pid(i, 1);  A[ p ] [ p ] = 1. ;  b[ p ] =  1.;
	p = pid(i,ny);  A[ p ] [ p ] = 1. ;  b[ p ] = -1.;
      }
    jLOOP
      {
	int p,i;
	p = pid(1, j);  A[ p ] [ p ] = 1. ;  b[ p ] =  1.;
	p = pid(nx,j);  A[ p ] [ p ] = 1. ;  b[ p ] = -1.;
      }
  }

  //  ==
  //  ||
  //  ||  Utility routines
  //  ||
  //  ==

  int pid(int i,int j) { return i + (j-1)*nx; }  // Given i-j, return point ID.

  #include "plotter.h"
  #include "gauss_seidel.h"

};



//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==

int main()
{

  cout << "\n";
  cout << "---------------------------------------------\n";
  cout << "|                                           |\n";
  cout << "| F I N I T E   D I F F D E R E N C E       |\n";
  cout << "| D E M O   C O D E                         |\n";
  cout << "|                                           |\n";
  cout << "---------------------------------------------\n";
  cout << "\n";

  LaplacianOnGrid F(10,10, "F" );

  cout << "                                          \n";
  cout << "Form Linear System:                       \n";
  cout << "------------------------------------------\n";
  cout << "                                          \n";

  F.FormLS();

  cout << "                                          \n";
  cout << "Solve Linear System:\n";
  cout << "------------------------------------------\n";
  cout << "                                          \n";

  F.GaussSeidel(200, F.b , F.phi );

  cout << "                                          \n";
  cout << "Plot Results:\n";
  cout << "------------------------------------------\n";
  cout << "                                          \n";

  F.plot("phi",F.phi);      

}
