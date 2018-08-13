// =====================================================================
// ||                                                                 ||
// ||                                                                 ||
// ||            S G H   S O L V E R                                  ||
// ||                                                                 ||
// ||   Created by: Scott R. Runnels, Ph.D.                           ||
// ||               Computational Physics Division                    ||
// ||               Los Alamos National Laboratoray                   ||
// ||                                                                 ||
// ||   For use in CVEN 5838, Computational Multi-Physics Production  ||
// ||   Software Development.                                         ||
// ||                                                                 ||
// =====================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include "stdio.h"
#include "math.h"
#include "./mat.h"
using std :: string;
using std :: vector;
using std :: stringstream;
using std :: cout;
using std :: endl;

#define pLOOP for ( int p = 0 ; p <= npts-1 ; ++p)
#define p_int_LOOP for ( int p = 1 ; p <= npts-2 ; ++p)
#define cLOOP for ( int c = 0 ; c <= ncells-1 ; ++c)
#define  D_  double*
#define _D_  double

int npts = 900;         int ncells = npts - 1;
_D_ dx = 0.;            _D_ length = 6.;  
_D_ dt = 1.e-03;        _D_ endtime = 10.;
_D_ density = 8.93;
_D_ energy = 0.;   
_D_ q2 = 0.1;
_D_ q1 = 0.1;
_D_ PistonSpeed = 0.01;

class State
{
public:

  Mat M; 
  _D_ p_x[1000]; _D_ p_v[1000]; _D_ c_x[1000];  _D_ c_a     [1000];
  _D_ c_r[1000]; _D_ c_E[1000]; _D_ c_p[1000];  _D_ c_ss    [1000];

  State(void)
  {  
    pLOOP p_x[p] = dx*float(p);       pLOOP p_v[p] = 0.;
    cLOOP c_r[c] = density;           cLOOP c_E[c] = energy;
    
    _D_ tmp_p = M.Pressure(c_r[0],c_E[0]);
                                          
    cLOOP c_p[c] = tmp_p ; cLOOP c_a[c]      = 0. ;
    cLOOP c_x[c] = 0.    ;
    cLOOP c_x[c] = (p_x[c]+p_x[c+1])/2.;           
  }                                                

  ~State(void){}

  void Update()
  {
    cLOOP c_p [c] = M.Pressure(c_r[c],c_E[c]);
    cLOOP c_ss[c] = M.SoundSpeed(c_r[c],c_E[c]);
  } 
};                                                                  


void SGH_Advance(State &Sold, D_ c_p, _D_ dt, State &Snew)
{                                                                  
  _D_ p_m = density * dx;                                          
  _D_ p_f[1000]; pLOOP p_f[p] = 0.;                                
                                                                   
  Snew.p_v[0] = PistonSpeed;                                       
  Snew.p_v[npts-1] = 0.;

  _D_ du, ss, rho, q, e_rate;
                                                                   
  cLOOP                                                            
    {                                                              
      du  = Sold.p_v[c+1] - Sold.p_v[c];                       
      rho = p_m/( Sold.p_x[c+1] - Sold.p_x[c] );              
      q   = 0.;                                                
      q  += q1 * rho * fabs(du) * Sold.c_ss[c];
      q  += q2 * rho * du            * du ;                         

      if (du > 0. )  q = 0.;                                       

      c_p[c  ] +=  q     ;                                                 
      p_f[c  ] += -c_p[c];                                         
      p_f[c+1] +=  c_p[c];                                         
    }                                                              
                                                                   
  p_int_LOOP Snew.p_v[p] =   Sold.p_v[p] + dt * p_f[p]/p_m;         
  pLOOP      Snew.p_x[p] =   Sold.p_x[p] + dt * Snew.p_v[p];             
  cLOOP      Snew.c_x[c] = ( Snew.p_x[c] + Snew.p_x[c+1] ) / 2.;        
                                                                   
  cLOOP                                                            
    {                                                              
      e_rate             = c_p[c]*Snew.p_v[c] - c_p[c]*Snew.p_v[c+1];
      Snew.c_E[c]        = Sold.c_E[c] +  dt * e_rate/p_m;           
      Snew.c_r[c]        = p_m / ( Snew.p_x[c+1] - Snew.p_x[c] );    
    }

  return;
}



int main() 
{          
  Mat M;   
           
  dx = length / float(npts - 1);

  _D_ time;
  State S; 
  State Smid;

  time = 0.;

  while ( time < endtime)
    {
      time += dt;

      SGH_Advance(S, S.c_p   , dt/2. , Smid);  Smid.Update(); 
      SGH_Advance(S, Smid.c_p, dt    , S   );  S.Update();

    }

  std::fstream f;
  f.open("e.plt",std::ios::out); 
  cLOOP  f << S.c_x[c] << " " << S.c_E[c] << endl;
  f.close(); 

}



