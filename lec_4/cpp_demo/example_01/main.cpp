#include "sayHello.h"

#define VD vector<double>

int main() 
{
  int err;
  err = sayHello();

  return err;
}


int sayHello()
{
  cout << "Hello" << endl;

  VD x_array;

  x_array.resize(10); // this array has 10 elements

  for(int i=0 ; i<10 ; ++i)
    {
      x_array[i] = i*i;
      cout << "x_array = " << x_array[i] << endl;
    }

  return 0;
}
