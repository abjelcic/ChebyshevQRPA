function Omegas = function_generateRPAfrequencies( n , Omegamax , maxdegblock )

    assert( length(n)           == 1           , 'n must be non-negative integer!'                                );
    assert( isreal(n)           == 1           , 'n must be non-negative integer!'                                );
    assert( floor(n)            == n           , 'n must be non-negative integer!'                                );
    assert( n                   >= 0           , 'n must be non-negative integer!'                                );
    assert( length(maxdegblock) == 1           , 'maxdegblock must be non-negative integer less or equal than n!' );
    assert( isreal(maxdegblock) == 1           , 'maxdegblock must be non-negative integer less or equal than n!' );
    assert( floor(maxdegblock)  == maxdegblock , 'maxdegblock must be non-negative integer less or equal than n!' );
    assert( maxdegblock         >= 0           , 'maxdegblock must be non-negative integer less or equal than n!' );
    assert( maxdegblock         <= n           , 'maxdegblock must be non-negative integer less or equal than n!' );
    assert( length(Omegamax)    == 1           , 'Omegamax must be non-negative real number!'                     );
    assert( isreal(Omegamax)    == 1           , 'Omegamax must be non-negative real number!'                     );
    assert( Omegamax            >= 0           , 'Omegamax must be non-negative real number!'                     );
    
    
    Omegas = [];
    while( length(Omegas) < n )
        nb = randi(maxdegblock);
        Omegas = [ Omegas , Omegamax*rand()*ones(1,nb) ];
    end
    Omegas = Omegas(1:n);
    
    return;
end
