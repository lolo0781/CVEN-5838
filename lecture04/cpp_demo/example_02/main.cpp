#include "sayHello.h"

#define VD vector<double>
#define iLOOP for(int i=0 ; i<10 ; ++i)


class Circle
{

public:

  // constructor
  Circle()
  {
    cout << "Hello, I am a circle that was just born" << endl;
  }

  // constructor that accepts args
  Circle(float init_radius, string _name)
  {
    r = init_radius;
    name = _name;
    cout << "Hello, I am a circle that was just born" << endl;
    cout << "  my initial radius is " << r << endl;
    cout << "  my name is " << name << endl;
  }

  // destructor
  ~Circle()
  {
    cout << "Goodbye, I am a circle that just died" << endl;
  }
  
  float xc, yc, r; // coords of center and radius
  string name;

  void printInfo()
  {
    cout << "The radius of my circle is " << r << endl;
  }

};


int main() 
{
  Circle C1, C2(34, "C2"); // Declare a variable of type Circle
  
  C1.xc = 0.;
  C1.yc = 0.;
  C1.r  = 1.;
  C1.printInfo();

  C2.xc = 0.;
  C2.yc = 0.;
  C2.printInfo();

  return 0;
}
