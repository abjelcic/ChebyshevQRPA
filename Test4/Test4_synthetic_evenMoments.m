clear;
clc;
close all;
addpath('../functions');
tic;

Np     = 1000;
Omegab = 250;


ndegenf = max([1,floor(Np*0.001)]);
Omegas1 = function_generateRPAfrequencies( floor(Np/2) , 200 , ndegenf );
Omegas2 = function_generateRPAfrequencies( ceil(Np/2)  , 50  , ndegenf );
Omegas  = [ Omegas1 , Omegas2 ];
iFO     = randn(1,Np) + 1j*randn(1,Np);
OFi     = conj(iFO);

[ A , B , F20 , F02 ] = function_generateRPAmatrices( Omegas , iFO , OFi , false );


PlotData = cell(0,1);

for kthMoment = [ 0 , 2 , 4 , 6 , 8 ]
    
    if( kthMoment == 0 )
        kernel = 'Jackson';
    else
        kernel = 'Dirichlet';
    end
    
    x = [];
    y = [];
    
    for N_it = arrayfun( @(n)2^n , [3:12] )
 
        mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , N_it );
        mun = function_applyKernel( mun , kernel );

        mk_Cheb    = function_ChebyshevMoment( mun , kthMoment , Omegab );

        mk_correct = function_kthMoment( Omegas , abs(iFO).^2 , kthMoment );

        relerr     = abs(mk_Cheb-mk_correct)/abs(mk_correct);

        fprintf('k = %1d, N_it = %5d, mk(Chebyshev) = %e, mk(correct) = %e, rel. err. = %e\n', kthMoment , N_it , mk_Cheb , mk_correct , relerr );

        x = [ x , N_it   ];
        y = [ y , relerr ];
    end
    
    PlotData{end+1} = { kthMoment , x , y };
    
    fprintf('\n');
end



figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
Legends  = cell(0);
AllMarks = {"o","d","^","s","x",".","+","*","v",">","<","p","h"};
for i = 1 : length(PlotData)
    kthMoment = PlotData{i}{1};
    N_it      = PlotData{i}{2};
    relerr    = PlotData{i}{3};
    
    Legends{i} = strcat("$k =$ ",num2str(kthMoment));
    
    semilogy( log2(N_it) , relerr , strcat(AllMarks{mod(i,length(AllMarks))},"-") , 'LineWidth' , 2 , 'MarkerSize' , 8 ); hold on;
    
    XtickLabels = cell(0);
    for i = 1 : length(N_it)
        %XtickLabels{i} = '2^{' + string( round( log2(N_it(i)) ) ) + '}';
        XtickLabels{i} = round( N_it(i) );
    end
    xticks( log2(N_it) );
    xticklabels( XtickLabels );
end
grid on; %grid minor;

legend(Legends,'Interpreter','latex');
xlabel('$N_{\mathrm{it}}$','Interpreter','latex');
ylabel('$|m_k-m_k^{\mathrm{(Chebyshev)}}|/m_k$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);    
    
time = toc;
fprintf( 'Total time: %.2f s.\n' , time );
