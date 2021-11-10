function dBdw = function_solveLinearResponseRPA( A , B , omega , gamma , F20 , F02 )

    assert( size(A,1)          == size(A,2)   , 'A must be a square matrix!'                                         );
    assert( size(B,1)          == size(B,2)   , 'B must be a square matrix!'                                         );
    assert( length(A)          == length(B)   , 'A and B must be of the same size!'                                  ); 
    assert( size(F20,2)        == 1           , 'F20 must be a one-column vector!'                                   );
    assert( size(F02,2)        == 1           , 'F02 must be a one-column vector!'                                   );
    assert( length(F20)        == length(F02) , 'F02 and F20 must be of the same size!'                              );
    assert( length(F20)        == length(A)   , 'F20 and F02 vectors must be of the same order as matrices A and B!' );
    assert( length(omega)      == 1           , 'omega must be real number!'                                         ); 
    assert( isreal(omega)      == 1           , 'omega must be real number!'                                         ); 
    assert( length(gamma)      == 1           , 'gamma must be real positive number!'                                ); 
    assert( isreal(gamma)      == 1           , 'gamma must be real positive number!'                                ); 
    assert( gamma              >  0           , 'gamma must be real positive number!'                                ); 
    
    n = length(F20);

    S = [ A      , B        ; conj(B)  , conj(A) ];
    X = [ eye(n) , zeros(n) ; zeros(n) , -eye(n) ];
    
    XY = ( S - (omega+1j*gamma)*X ) \ [-F20;-F02];
    
    Xw = XY(   1:n   );
    Yw = XY( n+1:n+n );
    
    Sw = [ F20 ; F02 ]' * [ Xw ; Yw ];
    
    dBdw = -1/pi * imag(Sw);
    
    return;
    
end

