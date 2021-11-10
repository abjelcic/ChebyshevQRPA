function Ink = function_Ink( n , k , Omega0 , Omegab )
    % Calculates 2/pi * integral_{[0,Omega0/Omegab]} x^k/sqrt(1-x^2) * Tn(x) dx
    
    assert( length(Omegab)     == 1      , 'Omegab must be real positive number!'    ); 
    assert( isreal(Omegab)     == 1      , 'Omegab must be real positive number!'    ); 
    assert( Omegab             >  0      , 'Omegab must be real positive number!'    ); 
    assert( length(Omega0)     == 1      , 'Omega0 must be real positive number!'    ); 
    assert( isreal(Omega0)     == 1      , 'Omega0 must be real positive number!'    ); 
    assert( Omega0             >  0      , 'Omega0 must be real positive number!'    );
    assert( Omega0             <= Omegab , 'Omega0 must be less or equal to Omegab!' );
    assert( length(n)          == 1      , 'n must be non-negative integer!'         );
    assert( isreal(n)          == 1      , 'n must be non-negative integer!'         );
    assert( n                  >= 0      , 'n must be non-negative integer!'         );
    assert( floor(n)           == n      , 'n must be non-negative integer!'         );
    assert( length(k)          == 1      , 'k must be non-negative integer!'         );
    assert( isreal(k)          == 1      , 'k must be non-negative integer!'         );
    assert( k                  >= 0      , 'k must be non-negative integer!'         );
    assert( floor(k)           == k      , 'k must be non-negative integer!'         );

    
    Ink = 0;
    for j = 0 : k
        if( mod(k-j,2) == 0 )
            fac1 = sinc(abs(n+j)/2);
            fac2 = sinc(abs(n-j)/2);
            fac3 = sinc( abs(n+j) * acos(Omega0/Omegab)/pi ) * acos(Omega0/Omegab)/pi;
            fac4 = sinc( abs(n-j) * acos(Omega0/Omegab)/pi ) * acos(Omega0/Omegab)/pi;
            fac5 = nchoosek( k , (k-j)/2 ) / 2^k;
            
            if( j == 0 )
                Ink = Ink + fac5 * 0.5 * ( fac1 + fac2 - 2*fac3 - 2*fac4 ); 
            else
                Ink = Ink + fac5 * 1.0 * ( fac1 + fac2 - 2*fac3 - 2*fac4 ); 
            end
        end
    end

    return;

end


