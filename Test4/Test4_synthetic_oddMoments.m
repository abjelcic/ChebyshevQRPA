clear;
clc;
close all;
addpath('../functions');

Np = 1000;


ndegenf = max([1,floor(Np*0.001)]);
Omegas1 = function_generateRPAfrequencies( floor(Np/2) , 200 , ndegenf );
Omegas2 = function_generateRPAfrequencies( ceil(Np/2)  , 50  , ndegenf );
Omegas  = [ Omegas1 , Omegas2 ];
iFO     = randn(1,Np) + 1j*randn(1,Np);
OFi     = conj(iFO);

[ A , B , F20 , F02 ] = function_generateRPAmatrices( Omegas , iFO , OFi , false );




for kthMoment = [ -9 , -7 , -5 , -3 , -1 , +1 , +3 , +5 , +7 , +9 ]

    S = [ A , B ; -conj(B) , -conj(A) ];
    
    xy = [ F20 ; -F02 ];
    for i = 1 : abs(kthMoment)
        if( kthMoment > 0 )
            xy = S * xy;
        else
            xy = S \ xy;
        end
    end
    mk_calculated = real( 0.5 * [ F20 ; F02 ]' * xy );  
    mk_correct    = function_kthMoment( Omegas , abs(iFO).^2 , kthMoment );
    
    relerr        = abs(mk_calculated-mk_correct)/abs(mk_correct);
    
    fprintf('k = %+2d, mk(calculated) = %e, mk(correct) = %e, rel. err. = %e\n' , kthMoment , mk_calculated , mk_correct , relerr );

end
