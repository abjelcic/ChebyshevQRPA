function [x,y] = function_LorentzianSmearing( Omegas , absiFO2 , absOFi2 , gamma_smear , Omegab , NevalPoints )
    
    assert( isvector(Omegas)    == 1               , 'Omegas must be a vector!'                     );
    assert( isreal(Omegas)      == 1               , 'Omegas must be a real vector!'                );
    assert( sum(Omegas>0)       == length(Omegas)  , 'Omegas must be positive vector!'              );
    assert( isvector(absiFO2)   == 1               , 'absiFO2 must be a vector!'                    );
    assert( isreal(absiFO2)     == 1               , 'absiFO2 must be a real vector!'               );
    assert( sum(absiFO2>0)      == length(absiFO2) , 'absiFO2 must be positive vector!'             );
    assert( length(Omegas)      == length(absiFO2) , 'Omegas an absiFO2 must be of the same length' );
    assert( isvector(absOFi2)   == 1               , 'absOFi2 must be a vector!'                    );
    assert( isreal(absOFi2)     == 1               , 'absOFi2 must be a real vector!'               );
    assert( sum(absOFi2>0)      == length(absOFi2) , 'absOFi2 must be positive vector!'             );
    assert( length(Omegas)      == length(absOFi2) , 'Omegas an absOFi2 must be of the same length' );
    assert( length(gamma_smear) == 1               , 'gamma_smear must be real positive number!'    ); 
    assert( isreal(gamma_smear) == 1               , 'gamma_smear must be real positive number!'    ); 
    assert( gamma_smear         >  0               , 'gamma_smear must be real positive number!'    );    
    assert( length(Omegab)      == 1               , 'Omegab must be real positive number!'         ); 
    assert( isreal(Omegab)      == 1               , 'Omegab must be real positive number!'         ); 
    assert( Omegab              >  0               , 'Omegab must be real positive number!'         );     
    assert( length(NevalPoints) == 1               , 'NevalPoints must be real natural number!'     );
    assert( isreal(NevalPoints) == 1               , 'NevalPoints must be real natural number!'     );
    assert( NevalPoints         >= 1               , 'NevalPoints must be real natural number!'     );

    
    x = linspace( -Omegab , +Omegab , NevalPoints );
    y = zeros( 1 , length(x) );
    for i = 1 : length(x)
        omega = x(i);
        
        dBdw = 0;
        for j = 1 : length(Omegas)
            dBdw = dBdw + absiFO2(j) * (gamma_smear/pi) / ( (omega-Omegas(j))^2 + gamma_smear^2 );
            dBdw = dBdw - absOFi2(j) * (gamma_smear/pi) / ( (omega+Omegas(j))^2 + gamma_smear^2 );
        end
        
        y(i) = dBdw;
    end
    
    return;

end

