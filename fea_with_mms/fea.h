using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include "stdio.h"
#include "math.h"
using std :: string;
using std :: vector;
using std :: stringstream;
using std :: cout;
using std :: endl;

// Custom short-cuts for this particular code    

#define _D_ double                 
#define _I_ int                    
#define VD  vector<double>         
#define VI  vector<double>         
#define VS  vector<string>         
#define VDD vector<vector<double> >
#define VII vector<vector<int> >   
#define VTD vector<TestData >      

// Mesh

#define eLOOP  for ( int e = 1 ; e <= nelem ; ++e )
#define iLOOP  for ( int i = 1 ; i <= npx   ; ++i )
#define jLOOP  for ( int j = 1 ; j <= npy   ; ++j )

// General

#define kLOOP4 for ( int k = 1 ; k <= 4     ; ++k )
#define mLOOP4 for ( int m = 1 ; m <= 4     ; ++m )
#define nLOOP4 for ( int n = 1 ; n <= 4     ; ++n )

#define kLOOP2 for ( int k = 1 ; k <= 2     ; ++k )
#define mLOOP2 for ( int m = 1 ; m <= 2     ; ++m )
#define nLOOP2 for ( int n = 1 ; n <= 2     ; ++n )

#define kLOOPgauss for ( int k = 1 ; k <= nGauss     ; ++k )

// Matrix

#define rowLOOP   for ( int row = 1 ; row <= nrows     ; ++row )
#define colLOOP   for ( int col = 1 ; col <= Alen[row] ; ++col )
#define colLOOPp1 for ( int col = 2 ; col <= Alen[row] ; ++col )
#define colLOOPbw for ( int col = 1 ; col <= bandwidth ; ++col )


