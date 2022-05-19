clear;
clc;
close all;
addpath('../functions');
tic;

fprintf('Readimg A.mat...\n'); A = function_readmatd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "A.mat" );
fprintf('Readimg B.mat...\n'); B = function_readmatd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "B.mat" );

[ Omegas , RPA_BIS , RPA_BIV ] = function_readSkyrmeRPAeigenfreq( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/" , "skyrme_rpa.out" );




fprintf('\nISOSCALAR\n');
F20 = function_readvecd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "F02_ISOSCALAR.vec" );
F02 = function_readvecd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "F20_ISOSCALAR.vec" );

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
    
    mk_correct    = function_kthMoment( Omegas , RPA_BIS , kthMoment );
    
    relerr        = abs(mk_calculated-mk_correct)/abs(mk_correct);
    
    fprintf('k = %+2d, mk(calculated) = %e, mk(correct) = %e, rel. err. = %e\n' , kthMoment , mk_calculated , mk_correct , relerr );

end




fprintf('\nISOVECTOR\n');
F20 = function_readvecd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "F02_ISOVECTOR.vec" );
F02 = function_readvecd( "../Test2/skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "F20_ISOVECTOR.vec" );

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
    
    mk_correct    = function_kthMoment( Omegas , RPA_BIV , kthMoment );
    
    relerr        = abs(mk_calculated-mk_correct)/abs(mk_correct);
    
    fprintf('k = %+2d, mk(calculated) = %e, mk(correct) = %e, rel. err. = %e\n' , kthMoment , mk_calculated , mk_correct , relerr );

end

time = toc;
fprintf( 'Total time: %.2f s.\n' , time );
