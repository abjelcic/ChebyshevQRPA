function mk = function_kthMoment( Omegas , absiF02 , k )

    assert( isvector(Omegas)    == 1               , 'Omegas must be a vector!'                     );
    assert( isreal(Omegas)      == 1               , 'Omegas must be a real vector!'                );
    assert( sum(Omegas>0)       == length(Omegas)  , 'Omegas must be positive vector!'              );
    assert( isvector(absiF02)   == 1               , 'absiF02 must be a vector!'                    );
    assert( isreal(absiF02)     == 1               , 'absiF02 must be a real vector!'               );
    assert( sum(absiF02>0)      == length(absiF02) , 'absiF02 must be positive vector!'             );
    assert( length(Omegas)      == length(absiF02) , 'Omegas an absiF02 must be of the same length' );   
    assert( length(k)           == 1               , 'k must be real integer number!'               );
    assert( isreal(k)           == 1               , 'k must be real integer number!'               );
    assert( floor(k)            == k               , 'k must be real integer number!'               );
    
    mk = 0;
    for i = 1 : length(Omegas)
        mk = mk + absiF02(i)*Omegas(i)^k;
    end
    
    return;
end

