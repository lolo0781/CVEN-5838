
#include "fea.h"


class LaplacianOnGrid
{
public:

  // Mesh
  // --------------------------------------------------------------------------
  int nnode, nelem;                    // No. nodes and elements
  VD  x,y;                             // Coordinates of each node
  _D_ length, dx, dy;                  // Size of domain and x/y spacing
  int npx, npy;                        // No. points   in x/y directions
  int nex, ney;                        // No. elements in x/y directions
  VII elem;                            // 2D array storing nodes in each element
                                       // n[e,pm] array from the text
  // A x = b
  // --------------------------------------------------------------------------
  VDD Acoef;  VD Soln;   VD b;         // The linear system in compressed format
  VII Jcoef; 
  VI  Alen;                            // No. of non-zeros in each row of A
  int bandwidth;                       // Max. no. of non-zeros in each row of A
  int nrows;                           // No. rows in the linear system

  // Element matrix
  // --------------------------------------------------------------------------

  VDD ematrix;                         // Element matrix storing all integrals on
                                       // the master element.
  int nGauss;                          // No. of Gauss quad. points.
  VDD GaussPts;                        // Gauss sampling points on the master element
  VDD GaussWts;                        // Gauss weights

  VD  dphi_dxi;                        // M.Element basis function derivaties (xi-eta plane)
  VD  dphi_deta;                       //
  VD  dw_dx;                           // Physical basis function derivaties (x-y plane)
  VD  dw_dy;                           //
  _D_ detJac;                          // Jacobian determinant magnitude

  int pid( int i , int j ) { return ( i + (j-1)*npx); }

  LaplacianOnGrid( int ncell_x  , int ncell_y )
  {
    // Mesh metadata
    // ------------------------------------------------------------------------

    cout << "(o) Mesh meta data\n";

    length = 1.;
    nex = ncell_x       ; ney = ncell_y      ;  nelem = nex*ney;
    dx  = length / nex  ;  dy = length / ney ;
    npx = nex + 1       ; npy = ney + 1      ;

    nnode     = npx * npy;
    bandwidth = 9;
    nrows     = nnode;

    // Mesh point locations
    // ------------------------------------------------------------------------

    cout << "(o) Mesh points\n";

    x.resize(nnode + 1)  ;  y.resize(nnode + 1);   // Physical point locations

    iLOOP jLOOP 
    {
      int p = pid(i,j);
      x[p]  = (i-1)*dx;
      y[p]  = (j-1)*dy;
    }

    // Which nodes (points) are in each element?
    // ------------------------------------------------------------------------

    cout << "(o) Mesh element formation \n";

    int ec = 0;   // Counts which element we are considering

    elem.resize(nelem + 1 );     eLOOP elem[e].resize( 4 + 1 );

    iLOOP jLOOP
    {
      if ( i < npx ) if ( j < npy )
      {
	++ec ; int p = pid(i,j);

	elem[ec][4] =  p + npx   ;   elem[ec][3] =  p + 1 + npx   ;
	elem[ec][1] =  p         ;   elem[ec][2] =  p + 1         ;

      }
    }
    
    // Memory for element matrix
    // ------------------------------------------------------------------------

    cout << "(o) Element \n";

    ematrix.resize( 4 + 1 );  kLOOP4 emtrix[k].resize(4 + 1);

    dphi_dxi.resize ( 4 + 1 );
    dphi_deta.resize( 4 + 1 );

    dw_dx.resize( 4 + 1 );
    dw_dy.resize( 4 + 1 );


    // Linear System
    // ------------------------------------------------------------------------

    cout << "(o) Linear System \n";

    Acoef.resize( nnode + 1 )  ; rowLOOP Acoef[row].resize( bandwidth + 1 );
    Jcoef.resize( nnode + 1 )  ; rowLOOP Jcoef[row].resize( bandwidth + 1 );
    Alen.resize ( nnode + 1 )  ;
    b.resize    ( nnode + 1 )  ;
    Soln.resize ( nnode + 1 )  ;

    rowLOOP colLOOPbw Acoef[row][col] = 0.  ;
    rowLOOP colLOOPbw Jcoef[row][col] = 0   ;
    rowLOOP           Jcoef[row][ 1 ] = row ;
    rowLOOP           Alen [row]      = 1   ;
    rowLOOP           b    [row]      = 0.  ;
    rowLOOP           Soln [row]      = 0.  ;

  }
  

}



int main()
{


}

