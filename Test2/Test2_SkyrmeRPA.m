clear;
clc;
close all;
addpath('../functions');
tic;

N_its         = [ 200 , 400 , 800 , 1600 , 3200 , 6400 ];
% Must be NPlotMesh > max(N_its)
%NPlotMesh    = 300000; % Recommended value with MATLAB
NPlotMesh     = 6500;   % Recommended value with Octave
kernel        = 'Lorentz';
lambdaLorentz = 1.5;
Omegab        = 250;
ISIV          = "ISOVECTOR";


fprintf('Reading A.mat...  \n'); A   = function_readmatd( "./skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "A.mat"                    );
fprintf('Reading B.mat...  \n'); B   = function_readmatd( "./skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , "B.mat"                    );
fprintf('Reading F20.vec...\n'); F20 = function_readvecd( "./skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , strcat("F02_",ISIV,".vec") );
fprintf('Reading F02.vec...\n'); F02 = function_readvecd( "./skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/" , strcat("F02_",ISIV,".vec") );

[gamma_smear,J] = function_readSkyrmeRPAsmearingJ( "./skyrme_rpa/skyrme_rpa_MODIFIED/" , "skyrme_rpa.in" );


[RPA_Omegas,RPA_BIS,RPA_BIV] = function_readSkyrmeRPAeigenfreq( "./skyrme_rpa/skyrme_rpa_MODIFIED/" , "skyrme_rpa.out" );
scaleIS = max(abs(RPA_BIS( RPA_Omegas>=0 & RPA_Omegas<=50 )));
scaleIV = max(abs(RPA_BIV( RPA_Omegas>=0 & RPA_Omegas<=50 )));



if( strcmp(ISIV,'ISOSCALAR') )
    [xPlotMesh1,yPlotMesh1] = function_readSkyrmeRPAstrength( "./skyrme_rpa/skyrme_rpa_MODIFIED/" , "Plot_Bel_IS.dat" );
    scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));
else
    [xPlotMesh1,yPlotMesh1] = function_readSkyrmeRPAstrength( "./skyrme_rpa/skyrme_rpa_MODIFIED/" , "Plot_Bel_IV.dat" );
    scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));
end




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
%     for i = 1 : length(RPA_Omegas)
%        x = RPA_Omegas(i);
%        if( strcmp(ISIV,'ISOSCALAR') )
%            y = RPA_BIS(i)/scaleIS * scale1;
%        else
%            y = RPA_BIV(i)/scaleIV * scale1;
%        end
%        line([+x,+x],[0,+y],'LineWidth',0.05,'Color',[0.5 0 0.5]); hold on;
%        line([-x,-x],[0,-y],'LineWidth',0.05,'Color',[0.5 0 0.5]); hold on;
%     end
%     grid on; grid minor;

    legend1 = strcat( "True response ($\gamma$ = " , num2str(gamma_smear) , " $\mathrm{MeV}$)" );
    legend2 = strcat( "KPM ($N_{\mathrm{it}} = $ " , num2str(N_it)        , ")"                );
    legend({legend1,legend2},'Interpreter','latex');
    
    xlim([0,50]); ylim([0,scale1*1.1]);
    xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex');
    ylabel(strcat("$dB(\omega)/d\omega$ $[\mathrm{fm^{",num2str(2*J),"}/MeV}]$"),'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20);
  
    fprintf('N_it = %5d, calculatin finished, Chebyshev smearing: %9.4f MeV.\n' , N_it , gamma_Cheb );
        
    pause(0.1);    
    
end
    
time = toc;
fprintf( 'Total time: %.2f s.\n' , time );
