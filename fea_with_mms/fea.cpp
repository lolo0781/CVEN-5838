//  ================================================================================
//  ||                                                                            ||
//  ||              fea                                                           ||
//  ||              -------------------------------------------                   ||
//  ||              F I N I T E   E L E M E N T                                   ||
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
//  ||          Copyright: 2018 Scott Runnels                                     ||
//  ||                                                                            ||
//  ================================================================================

//  ================================================================================
//  ||                                                                            ||
//  || Preliminaries (capabilities included via the C++ #include statement        ||
//  ||                                                                            ||
//  ================================================================================


#include "fea.h"

//  ==
//  ||
//  ||      C L A S S:   L A P L A C I A N O N G R I D 
//  ||
//  ||      Does the following:
//  ||
//  ||      (1) Forms a linear system Ax=b representing a 2D
//  ||          finite element approximation to Laplace's
//  ||          equation on a square using bilinear elements.
//  ||      (2) Performs Gauss-Seidel iteration on that linear
//  ||          system up to a specified number of iterations
//  ||          or a specified tolerance.
//  ||    
//  == 

class LaplacianOnGrid
{

public:

  // Mesh
  // ----------------------------------------------------------
  int nnode, nelem  ;         // No. points and elements
  VD  x, y          ;         // x-y coordinates of each point
  _D_ length, dx, dy;         // Size of domain & grid spacing
  int npx,npy       ;         // No. points   in x/y directions
  int nex,ney       ;         // No. elements in x/y directions
  VII elem          ;         // List of points in each element

  // A x = b
  // ----------------------------------------------------------
  VDD Acoef;  VD Soln;  VD b;  // The linear system
  VII Jcoef;
  VI  Alen;                    // Length of each row
  int bandwidth;               // Max. no. non-zero entries on a row
  int nrows;                   // Number of rows in the matrix
  VD Soln_anl;                 // Analytical solutionn (MMS)

  // Element matrix
  // ----------------------------------------------------------
  VDD ematrix  ;              // Element matrix
  int nGauss   ;              // Number of Gauss quad. points.
  VDD GaussPts ;              // Gauss quadrature points
  VD  GaussWts ;              // Gauss quadrature weights
  VD  phi      ;              // Basis function values (xi-eta plane) (MMS)
  VD  dphi_dxi ;              // Basis function derivatives (xi-eta plane) 
  VD  dphi_deta;              // 
  VD  dw_dx    ;              // Basis functions derivatives (x-y plane) 
  VD  dw_dy    ;              // 
  _D_ detJac   ;              // |Jacobian|
  _D_ xphys    ;              // Physical x location, xphys =
			      // xphys(xi, eta) (MMS)
  _D_ yphys    ;              // Physical y location, similarly (MMS)
  VD  RHS      ;              // Right hand side vector for bilinear
			      // element, stores source term phi_i
			      // (MMS)

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  LaplacianOnGrid(int ncell_x, int ncell_y)
  {
    // Mesh metadata
    // ----------------------------------------------------------

    cout << "(o) Mesh meta data \n";

    length = 1.;

    nex = ncell_x       ; ney = ncell_y;   nelem = nex*ney;
    dx  = length / nex  ;  dy = length / ney;
    npx = nex + 1       ; npy = ney + 1;

    nnode     = npx*npy;
    bandwidth = 9;
    nrows     = nnode;

    // Mesh point locations
    // ----------------------------------------------------------

    cout << "(o) Mesh points \n";

    x.resize(nnode + 1);
    y.resize(nnode + 1);

    iLOOP jLOOP 
    { 
      int n = pid(i,j) ; 
      x[n] = (i-1)*dx  ; 
      y[n] = (j-1)*dy  ;
    }

    // Points in each element
    // ----------------------------------------------------------

    cout << "(o) Mesh elements , nelem = " << nelem << endl;

    elem.resize(nelem + 1);      eLOOP elem[e].resize(4 + 1);

    int ec = 0;

    iLOOP jLOOP
    { 
      if ( i < npx ) if ( j < npy )
	{
	  ++ec; int p = pid(i,j);

	  elem[ec][4] = p + npx;  elem[ec][3] = p + 1 + npx;
	  elem[ec][1] = p      ;  elem[ec][2] = p + 1      ;
	}
    }

    // Element matrix - bilinear element (4 points)
    // ----------------------------------------------------------

    cout << "(o) Element \n";

    ematrix.resize( 4 + 1 );  
    kLOOP4 ematrix[k].resize( 4 + 1);
    
    RHS.resize( 4 + 1 );  // allocate memory for RHS and basis fn vals
			  // (MMS)

    dphi_dxi.resize ( 4 + 1 ); 
    dphi_deta.resize( 4 + 1 ); 

    dw_dx.resize  ( 4 + 1 ); 
    dw_dy.resize  ( 4 + 1 ); 

    // Gauss quadrature (4 integration points)
    // ----------------------------------------------------------

    cout << "(o) Quadrature \n";

    #include "gauss_quadrature.h"


    // Linear System
    // ----------------------------------------------------------

    cout << "(o) Linear system\n";

    Acoef.resize(nnode+1); rowLOOP Acoef[row].resize(bandwidth+1);
    Jcoef.resize(nnode+1); rowLOOP Jcoef[row].resize(bandwidth+1);
    Alen.resize (nnode+1);
    b.resize    (nnode+1);
    Soln.resize (nnode+1); 

    rowLOOP colLOOPbw Acoef[row][col] = 0. ;  // Zero all entries of matrix
    rowLOOP colLOOPbw Jcoef[row][col] = 0  ;  // Zero all col pointers  
    rowLOOP           Jcoef[row][ 1 ] = row;  // Put diagonal in first entry
    rowLOOP Alen           [row]      = 1  ;  // Each row has one entry
    rowLOOP b              [row]      = 0. ;
    rowLOOP Soln           [row]      = 0. ;

    // (MMS) Analytical Solution
    // ----------------------------------------------------------

    Soln_anl.resize(nnode + 1);     // Allocate memory for analytical
				    // solution (MMS)

    // set analytical solution (MMS)
    iLOOP jLOOP
    {
      int p = pid(i,j);
      Soln_anl[p] = MMS_Solution( x[p], y[p] );
    }

  }

  //  ==
  //  ||
  //  ||  Form Linear System Ax = b
  //  ||
  //  ==

  void FormLS()
  {
    // Loop over all elements, assembling the linear system

    cout << "(o) Form Linear System " << endl;

    eLOOP
      {
	// Form element matrix
	Form_ematrix(e);

	// Assemble element matrix into global
	Assemble_ematrix_into_global(e);
      }

    // Apply boundary conditions
    // ----------------------------------------------------------

    cout << "(o) Apply BCs " << endl;
    
    iLOOP
      {
	int p;
	p = pid(i,  1);  b[ p ] =  MMS_Solution(x[p],y[p]);    Acoef[p][1] = 1.  ;   Alen[p] = 1;
	p = pid(i,npy);  b[ p ] =  MMS_Solution(x[p],y[p]);    Acoef[p][1] = 1.  ;   Alen[p] = 1;
      }
    jLOOP
      {
	int p;
	p = pid(1  ,j);  b[ p ] =  MMS_Solution(x[p],y[p]);    Acoef[p][1] = 1.  ;   Alen[p] = 1;
	p = pid(npx,j);  b[ p ] =  MMS_Solution(x[p],y[p]);    Acoef[p][1] = 1.  ;   Alen[p] = 1;
      }

    return;

  }


  //  ==
  //  ||
  //  ||  Form element matrix
  //  ||
  //  ==

  void Form_ematrix( int eid )
  {
    mLOOP4 nLOOP4 ematrix[m][m] = 0;
    mLOOP4        RHS[m];             // init RHS to zero (MMS)

    kLOOPgauss
      {
	// points in xi-eta plane where we sample the integrand
	_D_ xi  = GaussPts[k][1]; 
	_D_ eta = GaussPts[k][2]; 

	// Gauss quadratic weight at point (xi, eta)
	_D_ wt  = GaussWts[k];

	// Evaluate our integrand at (xi, eta) for element eid
	BilinearElement(eid, xi, eta);

	mLOOP4
	  nLOOP4
	  ematrix[m][n] += detJac * wt *
	                   ( (dw_dx[m] * dw_dx[n]) + (dw_dy[m] * dw_dy[m]) );

	// Loop over basis fns for RHS (MMS)
	mLOOP4
	  RHS[m] -= wt * detJac * ( phi[m] * MMS_Source(xphys, yphys) );
      }
  }


  //  ==
  //  ||
  //  ||  Assemble Matrix
  //  ||
  //  ==

  void Assemble_ematrix_into_global(int eid)
  {
    bool Extend;  // flag tells if we need to extend the current row
		  // to accomodate the column to be added (we start
		  // with diagonal)
    _D_  row;     // row of matrix we are affecting
    int  col2add; // actual column number we need to add
    _D_  val2add; // value we need to add to that column

    // Loop over the entries in the 4x4 element matrix
    mLOOP4
      nLOOP4
    {
      // This is i and j for the D_(ij) entry in the text
      row     = elem[eid][m];
      col2add = elem[eid][n];

      // This is D_(ij). Need to put it into A
      val2add = ematrix[m][n];

      // Search for the col2add
      colLOOP
	{
	  if(Jcoef[row][col] == col2add)
	    {
	      Acoef[row][col] += ematrix[m][n];
	      Extend = false;
	    }
	}

      // we didn't find col2add in this row. Extend row
      if(Extend)
	{
	  Alen[row] += 1;
	  Jcoef[row][Alen[row]] = col2add;
	  Acoef[row][Alen[row]] = val2add;
	}
    }

    // add RHS to b (MMS)
    mLOOP4
      {
	row     = elem[eid][m];
	b[row] += RHS[m];
      }
  }


  //  ==
  //  ||
  //  ||  Bi-Linear Element Basis Functions
  //  ||
  //  ==

  // Helper functions: 1-D basis functions and derivatives
  // ------------------------------------------------------------

  _D_ Lag_1  ( _D_ xi ) { return ( -0.5 * (xi - 1. ) ); } 
  _D_ Lag_2  ( _D_ xi ) { return (  0.5 * (xi + 1. ) ); }
  _D_ Lag_1_d( _D_ xi ) { return ( -0.5              ); }
  _D_ Lag_2_d( _D_ xi ) { return (  0.5              ); }

  void BilinearElement( int eid , _D_ xi , _D_ eta )
  {

    // Basis function values at xi-eta (MMS)
    // ----------------------------------------------------------

    phi[1] = Lag_1(xi) * Lag_1(eta);
    phi[2] = Lag_2(xi) * Lag_1(eta);
    phi[3] = Lag_2(xi) * Lag_2(eta);
    phi[4] = Lag_1(xi) * Lag_2(eta);

    // Physical loc corresponding to the given xi,eta (MMS)
    // ----------------------------------------------------------

    xphys = 0.;
    yphys = 0.;

    kLOOP4
      {
	int p = elem[eid][k];   // global point id for this point in me
	xphys += x[p] * phi[k];
	yphys += y[p] * phi[k];
      }

    // Basis function derivatives in the xi-eta plane
    // ----------------------------------------------------------

    dphi_dxi[1] = Lag_1_d(xi)*Lag_1(eta);
    dphi_dxi[2] = Lag_2_d(xi)*Lag_1(eta);
    dphi_dxi[3] = Lag_2_d(xi)*Lag_2(eta);
    dphi_dxi[4] = Lag_1_d(xi)*Lag_2(eta);

    dphi_deta[1] = Lag_1(xi)*Lag_1_d(eta);
    dphi_deta[2] = Lag_2(xi)*Lag_1_d(eta);
    dphi_deta[3] = Lag_2(xi)*Lag_2_d(eta);
    dphi_deta[4] = Lag_1(xi)*Lag_2_d(eta);


    // Derivatives for the mapping from xi-eta to x-y
    // ----------------------------------------------------------

    _D_ dx_dxi = 0;        _D_ dx_deta = 0;
    _D_ dy_dxi = 0;        _D_ dy_deta = 0;

    kLOOP4
      {
	// global point id
	int point = elem[eid][k];
	
	dx_dxi  += x[point] * dphi_dxi[k];
	dx_deta += x[point] * dphi_deta[k];
	dx_dxi  += y[point] * dphi_dxi[k];
	dx_deta += y[point] * dphi_deta[k];
      }


    // Derivatives for the inverse mapping from xi-eta to x-y
    // ----------------------------------------------------------

    _D_ dxi_dx, dxi_dy, deta_dx, deta_dy;

    // determinate of Jacobian
    detJac = (dx_dxi * dy_deta) - (dx_deta * dy_dxi);

    dxi_dx  =  dy_deta / detJac;
    dxi_dy  = -dx_deta / detJac;
    deta_dx = -dy_dxi  / detJac;
    deta_dy =  dx_dxi  / detJac;


    // Goal: Basis fn derivs in x-y plane, via chain rule
    // ----------------------------------------------------------

    kLOOP4 dw_dx[k] = (dphi_dxi[k] * dxi_dx) + (dphi_deta[k]*deta_dx);
    kLOOP4 dw_dy[k] = (dphi_dxi[k] * dxi_dy) + (dphi_deta[k]*deta_dy);

  }
  //  ==
  //  ||
  //  ||  Manufactured Solutions (MMS)
  //  ||
  //  ==

  _D_ MMS_Solution(_D_ x, _D_ y)
  {
    return ( pow(x,3) + pow(y,3) );
  }

  _D_ MMS_Source(_D_ x, _D_ y)
  {
    return ( 6.*x + 6.*y );
  }

  //  ==
  //  ||
  //  ||  Utility routines
  //  ||
  //  ==

  int pid(int i , int j ) { return ( i + (j-1)*npx ); }

  #include "gauss_seidel.h"
  #include "plotter.h"

};


//  ==
//  ||
//  ||  Main Program
//  ||
//  ==

int main()
{
  cout << "\n";
  cout << "---------------------------------------------\n";
  cout << "| Finite Element Demonstration Code         |\n";
  cout << "---------------------------------------------\n";
  cout << "\n";

  cout << endl;
  cout << "--------------------------" << endl;
  cout << "Form FEA Solver "           << endl;
  cout << "--------------------------" << endl;
  cout << endl;

  LaplacianOnGrid MyFEA(50,50);

  cout << endl;
  cout << "--------------------------" << endl;
  cout << "Form Linear System "        << endl;
  cout << "--------------------------" << endl;
  cout << endl;

  MyFEA.FormLS();

  cout << endl;
  cout << "--------------------------" << endl;
  cout << "Solve Linear System " << endl;
  cout << "--------------------------" << endl;
  cout << endl;

  MyFEA.GaussSeidel(1000,0.0001);

  MyFEA.plot("Soln",MyFEA.Soln);                       

}
