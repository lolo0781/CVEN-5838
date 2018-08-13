//  ============================================================
//  |  Finite Difference code                                  | 
//  |  CVEN-5838                                               |
//  |  Louie J. Long                                           |
//  |  August 2018                                             |
//  ============================================================

#include "fd_llong.h"

class LaplacianOnGrid
{

  /*

    CLASS: LaplacianOnGrid
   
    Does the following:
    
    (1) Form a linear system Ax=b representing a 2D finite difference
        approx to Laplace eqn on a square.
  
    (2) Perfom Gauss-Seidel iteration on that linear system up to a
        specified number of iterations or a specified tolerance.

  */

public:

  _I_ nx, ny, nrows;
  _D_ length, dx;
  VDD A; VD phi; VD b;
  string Name;

  // ----------------------------------------
  //  *** Constructor: Initialize values ***
  // ----------------------------------------

  LaplacianOnGrid(int ncell_x, int ncell_y, string _Name)
  {
    nx     = ncell_x + 1;
    ny     = ncell_y + 1;
    nrows  = nx+ny;
    length = 1;
    dx     = length/(nx - 1);
    Name   = _Name;

    A.resize(nrows + 1); rLOOP A[r].resize(nrows + 1);
    b.resize(nrows + 1);
    phi.resize(nrows + 1);
  }


  // -----------------------------------
  //  *** Form Linear System Ax = b ***
  // -----------------------------------

  void FormLS()
  {
    // ** Initialize linear system **

    rLOOP cLOOP A[r][c] = 0.;
    rLOOP b[r] = 0.;
    rLOOP phi[r] = 0.;

    // ** Form matrix entries for the interior grid pts **

    _D_ dx2 = dx*dx;
    
    iLOOPi
      jLOOPi
    {
      int p = pid(i,j);
      A[p][p]    = -4./dx2;
      A[p][p+1]  = 1./dx2;
      A[p][p-1]  = 1./dx2;
      A[p][p+nx] = 1./dx2;
      A[p][p-nx] = 1./dx2;
    }

    // ** Apply boundary conditions **

    iLOOP
      {
	int p,j;
	p = pid(i,  1);  A[p][p] = 1. ; b[p] =  1.;
	p = pid(i, ny);  A[p][p] = 1. ; b[p] = -1.;
      }

    jLOOP
      {
	int p,i;
	p = pid(1,  j);  A[p][p] = 1. ; b[p] =  1.;
	p = pid(nx, j);  A[p][p] = 1. ; b[p] = -1.;
      }
  }

  // --------------------------
  //  *** Utility Routines ***
  // --------------------------

  int pid(int i, int j) 
  { 
    // given i,j return point ID
    return i + (j-1)*nx; 
  }

  #include "plotter.h"
  #include "gauss_seidel.h"

};


// ------------------------
//  **** Main Program ****
// ------------------------

int main()
{

  cout << "\n";
  cout << "-----------------------------------------\n";
  cout << " Finite Difference Code                  \n";
  cout << "-----------------------------------------\n";
  cout << "\n";

  LaplacianOnGrid LOG( 10, 10, "LOG" );

  cout << "\n";
  cout << "-----------------------------------------\n";
  cout << " Form Linear System                      \n";
  cout << "-----------------------------------------\n";
  cout << "\n";  

  LOG.FormLS();

  cout << "\n";
  cout << "-----------------------------------------\n";
  cout << " Solve Linear System                     \n";
  cout << "-----------------------------------------\n";
  cout << "\n";

  LOG.GaussSeidel( 200, LOG.b, LOG.phi );

  cout << "\n";
  cout << "-----------------------------------------\n";
  cout << " Plot Results                            \n";
  cout << "-----------------------------------------\n";
  cout << "\n";

  LOG.plot( "phi", LOG.phi );
}
