function [ A , B , F20 , F02 ] = function_generateRPAmatrices( Omegas , iFO , OFi , test )
    
    assert( isvector(Omegas) == 1              , 'Omegas must be a vector!'                );
    assert( isreal(Omegas)   == 1              , 'Omegas must be a real vector!'           );
    assert( sum(Omegas>0)    == length(Omegas) , 'Omegas must be positive vector!'         );
    assert( isvector(iFO)    == 1              , 'iFO must be a vector!'                   );
    assert( length(iFO)      == length(Omegas) , 'iFO and Omegas must have the same size!' );
    assert( isvector(OFi)    == 1              , 'OFi must be a vector!'                   );
    assert( length(OFi)      == length(Omegas) , 'OFi and Omegas must have the same size!' );
    
    n = length(Omegas);
    
    
    [D,~] = qr( randn(n)+1j*randn(n) );
    [C,~] = qr( randn(n)+1j*randn(n) );
    
    xi = [];
    yi = [];
    while( length(xi)<n )
        nb = randi( max([1,floor(n*0.005)]) );
        
        thetas = 3*rand(nb,1);
        
        xi = [ xi ; cosh(thetas) ];
        yi = [ yi ; sinh(thetas) ];
    end
    xi = xi(1:n);
    yi = yi(1:n);
    
    X =      D  * diag(xi) * C;
    Y = conj(D) * diag(yi) * C;
    
    
    
    
    A = +( X*diag(Omegas)*X' + conj(Y)*diag(Omegas)*transpose(Y) );
    B = -( X*diag(Omegas)*Y' + conj(Y)*diag(Omegas)*transpose(X) );
    
    A = 0.5 * ( A + A'           );
    B = 0.5 * ( B + transpose(B) );
    
    
    
    
    F20F02 = [ X , -conj(Y) ; -Y , conj(X) ] * [ reshape(iFO,[n,1]) ; reshape(OFi,[n,1]) ];
    F20 = F20F02(   1 :   n );
    F02 = F20F02( n+1 : n+n );
    
    
    
    if( test )
        fprintf('function_generateRPAmatrices test, the following values should be zero:\n');
        fprintf('-----------------------------------------------------------------------\n');
        
        M  = [ +eye(n) , zeros(n) ; zeros(n) , -eye(n) ];
        XY = [ X       , conj(Y)  ; Y        , conj(X) ];
        S  = [ A       , B        ; conj(B)  , conj(A) ];
        O  = blkdiag( diag(Omegas) , -diag(Omegas) );
        F  = [ F20 ; F02 ];
        
        [ ~ , NotPositiveDef ] = chol(S);

        disp( norm( A - A' ) );
        disp( norm( B - transp(B) ) );
        disp( norm( XY *M*XY' - M ) );
        disp( norm( XY'*M*XY  - M ) );
        disp( norm( S*XY - M*XY*O ) );
        disp( norm( [X',Y']*F - reshape(iFO,[n,1]) ) );
        disp( norm( [transp(Y),transp(X)]*F - reshape(OFi,[n,1]) ) );
        disp( NotPositiveDef );
        
        fprintf('-----------------------------------------------------------------------\n');
    end
    
    return;
end




