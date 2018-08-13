// ==
// ||
// || Error Norms (1): RMS
// ||
// ==

_D_ RMS_error()
{
_D_ sum = 0.;
int nval = 0;

iLOOPi jLOOPi
{
_D_ e = phi[ pid(i,j) ] - f_MMS( (i-1)*dx , (j-1)*dx );
sum += e*e;
++nval;
}
return sqrt(sum)/nval;
}

// ==
// ||
// || Error Norms (2): L2 error as an integral -- note, this one doesn't actually accomplish anything by subdividing each cell.
// ||
// ==

_D_ L2_error()
{
_D_ dA = (dx/10.) * (dx/10.);
_D_ area = 0.;
_D_ integral = 0.;

iLOOPi jLOOPi
{
for ( int K = 1 ; K <= 10 ; ++ K )
for ( int L = 1 ; L <= 10 ; ++ L )
{
_D_ e = phi[ pid(i,j) ] - f_MMS( (i-1)*dx , (j-1)*dx );
integral += e*e * dA;
area += 1. * dA;
}
}
return integral/area;
}

// ==
// ||
// || Error Norms (3): Max error 
// ||
// ==

_D_ max_error()
{
_D_ error = 0.;

iLOOPi jLOOPi
{
_D_ e = phi[ pid(i,j) ] - f_MMS( (i-1)*dx , (j-1)*dx );
e = fabs(e);
if ( e > error ) error = e;
}
return error;
}
