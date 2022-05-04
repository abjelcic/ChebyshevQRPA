clear;
clc;
close all;
addpath('../functions');

N_it          = 100000;
NPlotMesh     = 10000000;
kernel        = 'Lorentz';
lambdaLorentz = 1.5;



[xPlotMesh1,yPlotMesh1,gamma_smear] = function_readDIRQFAMstrength( "./DIRQFAM/DIRQFAMCheb/output/QFAM_output/" , "strength_Zr100_ISJ2K1.out" );
scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));


[mun,Omegab] = function_readDIRQFAMmun( "./DIRQFAM/DIRQFAMCheb/output/QFAM_output/" , "mu_Zr100_ISJ2K1.out" );
assert( length(mun) >= 2*N_it+1 , 'Not enough mun coefficients.' );
mun = mun( 1 : 2*N_it+1 );
mun = function_applyKernel( mun , kernel , lambdaLorentz );

gamma_Cheb = Omegab*lambdaLorentz/(2*N_it+1);

[xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( mun , Omegab , NPlotMesh );
scale2 = max(abs(yPlotMesh2( xPlotMesh2>=0 & xPlotMesh2<=50 )));



figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
plot( xPlotMesh1 , yPlotMesh1 , 'r.-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
plot( xPlotMesh2 , yPlotMesh2 , 'b:'  , 'LineWidth' , 2.5 , 'MarkerSize' , 3 ); hold on;
%grid on; grid minor;

legend1 = strcat("True response ($\gamma$ = " , num2str(gamma_smear) , " $\mathrm{MeV}$)" );
legend2 = strcat("KPM ($N_{\mathrm{it}} = $ " , num2str(N_it)        , ")"                );
legend({legend1,legend2},'Interpreter','latex');

xlim([0,17]); ylim([0,+Inf]);
xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex');
ylabel('$dB(\omega)/d\omega$ $[\mathrm{fm^4/MeV}]$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
    
fprintf('N_it = %5d, Chebyshev smearing: %9.4f MeV.\n' , N_it , gamma_Cheb );

