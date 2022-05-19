clear;
clc;
close all;
addpath('../functions');
tic;

N_its         = [ 200 , 400 , 800 , 1600 , 3200 , 6400 ];
Np            = 1000;
% Must be NPlotMesh > max(N_its)
%NPlotMesh    = 300000; % Recommended value with MATLAB
NPlotMesh     = 6500;   % Recommended value with Octave
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

fprintf('Generating RPA matrices.\n');
[ A , B , F20 , F02 ] = function_generateRPAmatrices( Omegas , iFO , OFi , false );








figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
for i = 1 : length(Omegas)
   x = +Omegas(i);
   y = +abs(iFO(i))^2;
   line([x,x],[0,y],'LineWidth',0.05,'Color',[0.5 0 0.5]); hold on;
   
   x = -Omegas(i);
   y = -abs(OFi(i))^2;
   line([x,x],[0,y],'LineWidth',0.05,'Color',[0.5 0 0.5]); hold on;
end
grid on; %grid minor;
% xlim([-Omegab,+Omegab]);
% ylim( [ -1.1*max([ abs(iFO).^2 , abs(OFi).^2 ]) , +1.1*max([ abs(iFO).^2 , abs(OFi).^2 ]) ] );
xlim([0,+Omegab]);
ylim( [ 0 , +1.1*max([ abs(iFO).^2 , abs(OFi).^2 ]) ] );
xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex');
ylabel('$dB(\omega)/d\omega$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
pause(0.01);  

fprintf('Calculating true response funtion with Lorentzian smearing.\n');
[xPlotMesh1,yPlotMesh1] = function_LorentzianSmearing( Omegas , abs(iFO).^2 , abs(OFi).^2 , gamma_smear , Omegab , NPlotMesh );
scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));




for N_it = N_its
    
    mun = function_ChebyshevCoefficients( A , B , F20 , F02 , Omegab , N_it );
    mun = function_applyKernel( mun , kernel , lambdaLorentz );
    
    gamma_Cheb = Omegab*lambdaLorentz/length(mun);
    
    fprintf( 'N_it = %5d, calculating Chebyshev coefficients.\n' , N_it);
    [xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( N_it , mun , Omegab , NPlotMesh );
    scale2 = max(abs(yPlotMesh2( xPlotMesh2>=0 & xPlotMesh2<=50 )));
    

    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,0.65,0.8]);
    plot( xPlotMesh1 , yPlotMesh1 , 'r.-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
    plot( xPlotMesh2 , yPlotMesh2 , 'b:'  , 'LineWidth' , 2.5 , 'MarkerSize' , 3 ); hold on;
    %grid on; grid minor;
    
    legend1 = strcat( "True response ($\gamma$ = " , num2str(gamma_smear) , " $\mathrm{MeV}$)" );
    legend2 = strcat( "KPM ($N_{\mathrm{it}} = $ " , num2str(N_it)        , ")"                );
    legend({legend1,legend2},'Interpreter','latex');
    
    xlim([0,50]); ylim([0,scale1*1.1]);
    xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex');
    ylabel('$dB(\omega)/d\omega$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20);
    
    fprintf('N_it = %5d, calculatin finished, Chebyshev smearing: %9.4f MeV.\n' , N_it , gamma_Cheb );
        
    pause(0.1);    
    
end

time = toc;
fprintf( 'Total time: %.2f s.\n' , time );
