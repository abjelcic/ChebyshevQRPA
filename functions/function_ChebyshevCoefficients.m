function mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , Nit )
    
    assert( size(A,1)          == size(A,2)   , 'A must be a square matrix!'                                         );
    assert( size(B,1)          == size(B,2)   , 'B must be a square matrix!'                                         );
    assert( length(A)          == length(B)   , 'A and B must be of the same size!'                                  ); 
    assert( size(F20,2)        == 1           , 'F20 must be a one-column vector!'                                   );
    assert( size(F02,2)        == 1           , 'F02 must be a one-column vector!'                                   );
    assert( length(F20)        == length(F02) , 'F02 and F20 must be of the same size!'                              );
    assert( length(F20)        == length(A)   , 'F20 and F02 vectors must be of the same order as matrices A and B!' );
    assert( length(Omegab)     == 1           , 'Omegab must be real positive number!'                               ); 
    assert( isreal(Omegab)     == 1           , 'Omegab must be real positive number!'                               ); 
    assert( Omegab             >  0           , 'Omegab must be real positive number!'                               ); 
    assert( length(Nit)        == 1           , 'Nit must be real natural number!'                                   );
    assert( isreal(Nit)        == 1           , 'Nit must be real natural number!'                                   );
    assert( Nit                >= 1           , 'Nit must be real natural number!'                                   );

    n = length(A);
    
    S = [ A , B ; conj(B) , conj(A) ];
    X = [ speye(n) , sparse(n,n) ; sparse(n,n) , -speye(n) ]; 
    F = [ F20 ; F02 ];
    
    alpha_old = X*F;
    alpha_new = 1/Omegab * (X*(S*(X*F)));
    
    mu0 = 0.5 * real( F' * alpha_old );
    mu1 = real( F' * alpha_new );
    
    mun = zeros( 2*Nit-1 , 1 );
    mun(0 +1) = mu0;
    for n = 1 : Nit
        % Now alpha1 contains alpha_{n}, and alpha0 contains alpha_{n-1} for given n
        mun(2*n-1 +1) = -  mu1 + 2 * real( alpha_old' * (X*alpha_new) );
        mun(2*n   +1) = -2*mu0 + 2 * real( alpha_new' * (X*alpha_new) );
        
        % Update alpha0 and alpha1
        if( n ~= Nit )
            alpha_tmp = alpha_new;
            alpha_new = 2/Omegab * (X*(S*alpha_new)) - alpha_old;
            alpha_old = alpha_tmp;
        end
    end
    
               
end





