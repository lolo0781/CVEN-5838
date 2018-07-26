    nGauss = 9;

    GaussPts.resize(nGauss + 1);  kLOOPgauss GaussPts[k].resize( 2 + 1 );
    GaussWts.resize(nGauss + 1); 

    _D_ xgA=  sqrt(3./5.)       ;   _D_ wgA = 5./9.;
    _D_ xgB=  0.                ;   _D_ wgB = 8./9.;

    //      xi                     eta                weight
    //  - - - - - - - -      - - - - - - - -       - - - - - - -

    GaussPts[1][1] = -xgA ; GaussPts[1][2] = -xgA ; GaussWts[1] = wgA;
    GaussPts[2][1] =  xgB ; GaussPts[2][2] = -xgA ; GaussWts[2] = wgA;
    GaussPts[3][1] =  xgA ; GaussPts[3][2] = -xgA ; GaussWts[3] = wgA;
    GaussPts[4][1] = -xgA ; GaussPts[4][2] =  xgB ; GaussWts[4] = wgA;
    GaussPts[5][1] =  xgB ; GaussPts[5][2] =  xgB ; GaussWts[5] = wgB;
    GaussPts[6][1] =  xgA ; GaussPts[6][2] =  xgB ; GaussWts[6] = wgA;
    GaussPts[7][1] = -xgA ; GaussPts[7][2] =  xgA ; GaussWts[7] = wgA;
    GaussPts[8][1] =  xgB ; GaussPts[8][2] =  xgA ; GaussWts[8] = wgA;
    GaussPts[9][1] =  xgA ; GaussPts[9][2] =  xgA ; GaussWts[9] = wgA;
