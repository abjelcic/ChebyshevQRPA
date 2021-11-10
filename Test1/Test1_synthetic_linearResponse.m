% MATLAB R2018a
clear;
clc;
close all;
addpath('../functions');

Np            = 50;
NPlotMesh     = 300000;
gamma_smear   = 0.05;
kernel        = 'Lorentz';
lambdaLorentz = 1.5;
Omegab        = 250;


ndegenf = max([1,floor(Np*0.001)]);
Omegas1 = function_generateRPAfrequencies( floor(Np/2) , 200 , ndegenf );
Omegas2 = function_generateRPAfrequencies( ceil(Np/2)  , 50  , ndegenf );
Omegas  = [ Omegas1 , Omegas2 ];
iFO     = randn(1,Np) + 1j*randn(1,Np);
OFi     = conj(iFO);

[ A , B , F20 , F02 ] = function_generateRPAmatrices( Omegas , iFO , OFi , false );









xPlotMesh1 = linspace( -Omegab , +Omegab , NPlotMesh );
yPlotMesh1 = zeros(1,NPlotMesh);
for i = 1 : NPlotMesh
    omega = xPlotMesh1(i);
    dBdw  = function_solveLinearResponseRPA( A , B , omega , gamma_smear , F20 , F02 );
    
    yPlotMesh1(i) = dBdw;
end
scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));





for N_it = [ 200 , 400 , 800 , 1600 , 3200 , 6400 ]
    
    mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , N_it );
    mun = function_applyKernel( mun , kernel , lambdaLorentz );
    
    gamma_Cheb = Omegab*lambdaLorentz/length(mun);
    
    [xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( mun , Omegab , NPlotMesh );
    scale2 = max(abs(yPlotMesh2( xPlotMesh2>=0 & xPlotMesh2<=50 )));

    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
    plot( xPlotMesh1 , yPlotMesh1 , 'r.-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
    plot( xPlotMesh2 , yPlotMesh2 , 'b:'  , 'LineWidth' , 2.5 , 'MarkerSize' , 3 ); hold on;
    %grid on; grid minor;
    
    legend1 = 'True response ($\gamma$ = ' + string(gamma_smear) + ' $\mathrm{MeV}$)';
    legend2 = 'KPM ($N_{\mathrm{it}} = $ ' + string(N_it) + ')';
    legend({legend1,legend2},'Interpreter','latex','FontSize',20);
    
    xlim([0,50]); ylim([0,scale1*1.1]);
    xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex','FontSize',30);
    ylabel('$dB(\omega)/d\omega$','Interpreter','latex','FontSize',30);
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20);
    
    fprintf('N_it = %5d, Chebyshev smearing: %9.4f MeV.\n' , N_it , gamma_Cheb );
        
    pause(0.1);    
    
end

