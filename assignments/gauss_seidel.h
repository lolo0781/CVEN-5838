  //  ==
  //  ||
  //  ||  Perform Gauss-Seidel iterations on nrows X nrows linear
  //  ||  system:  A*Solution = RHS
  //  ||
  //  ==

  void GaussSeidel(_I_ max_iter , VD RHS, VD &Solution)
  {
    _I_ converged, it_converged;
    _I_ iter = 0;
    _D_ newval;
    _D_ cur_delta = 0.;
    _D_ max_delta = 0.;
    _D_ tol       = 1.e-04;
    
    converged = 0;
    
    while ( converged == 0 )
      {
	// Check for excessive iterations

	if ( ++iter > max_iter )
	  { 
	    cout << "Gauss-Seidel in " << Name << ": max_iter of " << max_iter << " reached.\n";
	    return;
	  }

	// Initialize flags and values for convergence check
	
	max_delta    = 0.;
	it_converged = 1;

	// Perform one iteration of GS, loop over matrix rows

	rLOOP
	  {
	    // Compute new guess for row r

	    newval = RHS[r];
	    cLOOP if ( c != r ) newval -=  A[r][c] * Solution[c];
	    newval    /= A[r][r];

	    // Convergence check

	    cur_delta  = fabs(Solution[r] - newval);
	    if ( cur_delta > tol ) it_converged = 0;

	    // Record new value in solution

	    Solution[r]       = newval;
	  }

	// Make note of the convergence state
	
	converged = it_converged;
	
      }

    cout << "Gauss-Seidel in " << Name << ": converged in " << iter << " iterations.\n";
  }

