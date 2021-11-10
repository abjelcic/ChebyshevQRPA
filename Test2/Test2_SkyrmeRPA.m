% MATLAB R2018a
clear;
clc;
close all;
addpath('../functions');

NPlotMesh     = 300000;
kernel        = 'Lorentz';
lambdaLorentz = 1.5;
Omegab        = 250;
ISIV          = "ISOVECTOR";


A   = function_readmatd( './skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'A.mat'            );
B   = function_readmatd( './skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'B.mat'            );
F20 = function_readvecd( './skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F02_'+ISIV+'.vec' );
F02 = function_readvecd( './skyrme_rpa/skyrme_rpa_MODIFIED/out_matvec/' , 'F02_'+ISIV+'.vec' );

[gamma_smear,J] = function_readSkyrmeRPAsmearingJ( './skyrme_rpa/skyrme_rpa_MODIFIED/' , 'skyrme_rpa.in' );


[RPA_Omegas,RPA_BIS,RPA_BIV] = function_readSkyrmeRPAeigenfreq( './skyrme_rpa/skyrme_rpa_MODIFIED/' , 'skyrme_rpa.out' );
scaleIS = max(abs(RPA_BIS( RPA_Omegas>=0 & RPA_Omegas<=50 )));
scaleIV = max(abs(RPA_BIV( RPA_Omegas>=0 & RPA_Omegas<=50 )));



if( strcmp(ISIV,'ISOSCALAR') )
    [xPlotMesh1,yPlotMesh1] = function_readSkyrmeRPAstrength( './skyrme_rpa/skyrme_rpa_MODIFIED/' , 'Plot_Bel_IS.dat' );
    scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));
else
    [xPlotMesh1,yPlotMesh1] = function_readSkyrmeRPAstrength( './skyrme_rpa/skyrme_rpa_MODIFIED/' , 'Plot_Bel_IV.dat' );
    scale1 = max(abs(yPlotMesh1( xPlotMesh1>=0 & xPlotMesh1<=50 )));
end




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

    legend1 = 'True response ($\gamma$ = ' + string(gamma_smear) + ' $\mathrm{MeV}$)';
    legend2 = 'KPM ($N_{\mathrm{it}} = $ ' + string(N_it) + ')';
    legend({legend1,legend2},'Interpreter','latex','FontSize',20);
    
    xlim([0,50]); ylim([0,scale1*1.1]);
    xlabel('$\omega$ $[\mathrm{MeV}]$','Interpreter','latex','FontSize',30);
    ylabel('$dB(\omega)/d\omega$ $[\mathrm{fm^{'+string(2*J)+'}/MeV}]$','Interpreter','latex','FontSize',30);
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20);
  
    fprintf('N_it = %5d, Chebyshev smearing: %9.4f MeV.\n' , N_it , gamma_Cheb );
        
    pause(0.1);    
    
end
    
    