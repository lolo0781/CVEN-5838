#include "sayHello.h"

#define VD    vector<double>
#define VDD   vector< vector<double> >
#define iLOOP for(int i=0 ; i<10 ; ++i)


class Grid
{

  // What we want:

  //                                   (nx,ny)
  //        +----+----+----+----+----+----+
  //        |    |    |    |    |    |    |
  //        +----+----+----+----+----+----+
  //        |    |    |    |    |    |    |
  // (1,3)  +----+----+----+----+----+----+
  //        |    |    |    |    |    |    |
  // (1,2)  +----+----+----+----+----+----+
  //        |    |    |    |    |    |    |
  //        +----+----+----+----+----+----+
  //      (1,1) (2,1) ...               (nx,1)
  

public:

  Grid(int _nx, int _ny, double _x1, double _x2, double _y1, double _y2)
  {
    nx = _nx; ny = _ny;
    x1 = _x1; x2 = _x2;
    y1 = _y1; y2 = _y2;

    dx = (x2-x1)/nx;
    dy = (y2-y1)/ny;

    // allocate memory for vectors
    px.resize(nx); py.resize(ny);
    for(int i=0 ; i<nx ; ++i)
      {
	px[i].resize(ny);
	py[i].resize(nx);
      }

    for(int i=0 ; i<nx ; ++i)
      {
	for(int j=0 ; j<ny ; ++j)
	  {
	    px[i][j] = i*dx;
	    py[i][j] = j*dy;
	  }
      }
  }

  int nx, ny;
  double x1, x2, y1, y2;
  double dx, dy;

  VDD px;
  VDD py;

};


int main() 
{
  Grid G1(10, 10, 0., 1., 0., 1.);

  return 0;
}
