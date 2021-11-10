function [xk,fk] = function_fftEvaluateChebSeries( cn , Omegab , NevalPoints )
    % cn contains [ c_0 , c_1 , c_2 , ... , c_{N-1} ], i.e. cn(n +1) = c_n, for n = 0,1,2,...,N-1
    % Evaluates the Chebyshev series f(x) = (2/pi)/sqrt(Omegab^2-x^2) * sum_{n=0}^{N-1} c_n * T_n(x/Omegab) on
    % NevalPoints points in [-Omegab,+Omegab] interval.
    % Returns xk and fk, where xk(k +1) = Omegab * cos(pi*(k+1/2)/NevalPoints), and
    % fk(k +1) = f(xk(k +1)), for k=0,1,...,NevalPoints-1.
    
    assert( isvector(cn)        == 1            , 'cn must be a vector!'                      );
    assert( isreal(cn)          == 1            , 'cn must be a real vector!'                 );
    assert( length(Omegab)      == 1            , 'Omegab must be real positive number!'      ); 
    assert( isreal(Omegab)      == 1            , 'Omegab must be real positive number!'      ); 
    assert( Omegab              >  0            , 'Omegab must be real positive number!'      ); 
    assert( length(NevalPoints) == 1            , 'NevalPoints must be real natural number!'  );
    assert( isreal(NevalPoints) == 1            , 'NevalPoints must be real natural number!'  );
    assert( NevalPoints         >= 1            , 'NevalPoints must be real natural number!'  );
    assert( NevalPoints         >  2*length(cn) , 'NevalPoints must be at least 2xlength(cn)' );
    
    
    
    
    N = length(cn);
    Cn = zeros( 2*NevalPoints , 1 );
    for n = 1 : N
        Cn(n) = cn(n) * exp( -1j*pi*(n-1)/(2*NevalPoints) );
    end
    
    fk = real( fft(Cn) );
    fk = fk( (0)+1 : (NevalPoints-1)+1 );
    
    xk = Omegab .* cos( (pi/NevalPoints) .* ( [ 0 : NevalPoints-1 ]' + 0.5 ) );
    
    xk = reshape(xk,[1,length(xk)]);
    fk = reshape(fk,[1,length(fk)]);
    
    xk = fliplr(xk);
    fk = fliplr(fk);
    
    fk = 2/pi * 1 ./ sqrt( Omegab^2 - xk.^2 ) .* fk;
    
    return;
end

