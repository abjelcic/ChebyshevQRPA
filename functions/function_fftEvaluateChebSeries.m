function [wk,fk] = function_fftEvaluateChebSeries( N_it , mun , Omegab , NevalPoints )
    % mun contains [ mu_0 , mu_1 , mu_2 , ... , mu_{2*N_it} ], i.e. mun(n+1) = mu_n, for n = 0,1,2,...,2*N_it
    % Evaluates the Chebyshev series f(w) = (2/pi)/sqrt(Omegab^2-w^2) * sum_{n=0}^{N_it} mu_n * T_n(w/Omegab) on
    % NevalPoints points in [-Omegab,+Omegab] interval.
    % Returns wk and fk, where wk(k) = Omegab * cos(pi*(k-1/2)/NevalPoints), and
    % fk(k) = f(wk(k)), for k=1,...,NevalPoints.
    
    assert( isvector(mun)       == 1    , 'mun must be a vector!'                    );
    assert( isreal(mun)         == 1    , 'mun must be a real vector!'               );
    assert( length(Omegab)      == 1    , 'Omegab must be real positive number!'     ); 
    assert( isreal(Omegab)      == 1    , 'Omegab must be real positive number!'     ); 
    assert( Omegab              >  0    , 'Omegab must be real positive number!'     ); 
    assert( length(NevalPoints) == 1    , 'NevalPoints must be real natural number!' );
    assert( isreal(NevalPoints) == 1    , 'NevalPoints must be real natural number!' );
    assert( NevalPoints         >= 1    , 'NevalPoints must be real natural number!' );
    assert( NevalPoints         >  N_it , 'NevalPoints must be larget than N_it!'    );
    
    
    Mun = zeros( 2*NevalPoints , 1 );
    for n = 1 : 2*N_it+1
        Mun(n) = mun(n) * exp( -1j*pi * (n-1)/(2*NevalPoints) );
    end 
    
    fk = real( fft(Mun) );
    fk = fk( 1 : NevalPoints );
    
    wk = Omegab .* cos( (pi/NevalPoints) .* ( [1:NevalPoints]' - 0.5  ) );
    
       
    wk = reshape(wk,[1,length(wk)]);
    fk = reshape(fk,[1,length(fk)]);
    
    wk = fliplr(wk);
    fk = fliplr(fk);
    
    fk = (2/pi) ./ sqrt( Omegab^2 - wk.^2 ) .* fk;
    
    return;
end

