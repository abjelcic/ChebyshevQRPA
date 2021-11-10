function mk = function_ChebyshevMoment( mun , k , Omegab )
    
    assert( isvector(mun)      == 1      , 'mun must be real vector!'                );
    assert( isreal(mun)        == 1      , 'mun must be real vector!'                );
    assert( length(k)          == 1      , 'k must be non-negative integer!'         );
    assert( isreal(k)          == 1      , 'k must be non-negative integer!'         );
    assert( k                  >= 0      , 'k must be non-negative integer!'         );
    assert( floor(k)           == k      , 'k must be non-negative integer!'         );
    assert( length(Omegab)     == 1      , 'Omegab must be real positive number!'    ); 
    assert( isreal(Omegab)     == 1      , 'Omegab must be real positive number!'    ); 
    assert( Omegab             >  0      , 'Omegab must be real positive number!'    ); 
    
    
    Ink = zeros( length(mun) , 1 );
    for n = 0 : length(mun)-1
        Ink(n +1) = function_Ink( n , k , Omegab , Omegab );
    end
    
    mk = Omegab^k * dot( mun , Ink );

    return;

end


