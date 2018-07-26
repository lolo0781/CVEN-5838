  //  ==
  //  ||
  //  ||  Perform Gauss-Seidel iterations on nrows X nrows linear
  //  ||  system:  A*Solution = RHS
  //  ||
  //  ==

  void GaussSeidel(_I_ max_iter , _D_ tol )
  {
    _I_ converged, it_converged;
    _I_ iter = 0;
    _D_ newval;
    _D_ delta;
    
    converged = 0;
    
    while ( converged == 0 )
      {  
	// Iteration check

	if ( ++iter > max_iter ) {  cout << "Failed to converge\n";  return; }

	// Update each row

	it_converged = 1;

	rowLOOP    
	  {
	    newval            = b[row]; 
	    colLOOPp1 newval -= Acoef[row][col]*Soln[Jcoef[row][col]]; 
	    newval           /= Acoef[row][1];
	    delta             = fabs(Soln[row] - newval); 
	    Soln[row]         = newval;                                                                   
	    
	    if ( delta > tol ) it_converged = 0;
	  }
	
	converged = it_converged; 
      }

    cout << "Gauss-Seidel converged in " << iter << " iterations.\n";
    return;
  }

