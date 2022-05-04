clear;
clc;
close all;
addpath('../functions');


figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,1.0,1.0]);


    
N = 64;
mun = cos( pi/2 * [0:N-1] ) ./ ( 1 + ([0:N-1]==0) );

munJackson   = function_applyKernel( mun , 'Jackson'   );
munDirichlet = function_applyKernel( mun , 'Dirichlet' );


sigma = pi/N;
xPlotMesh1 = linspace(-1,+1,300);
yPlotMesh1 = 1/sqrt(2*pi*sigma^2) * exp( -0.5 * ( xPlotMesh1./sigma ).^2 );

[xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( munJackson   , 1.0 , 100000 );
[xPlotMesh3,yPlotMesh3] = function_fftEvaluateChebSeries( munDirichlet , 1.0 , 100000 );


plot( xPlotMesh1 , yPlotMesh1 , 'ro'  , 'LineWidth' , 2.0 , 'MarkerSize' , 3 ); hold on;
plot( xPlotMesh2 , yPlotMesh2 , 'b-'  , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
plot( xPlotMesh3 , yPlotMesh3 , 'k--' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;


legend1 = "Gaussian $\sigma = \pi/N$";
legend2 = strcat( "KPM ($N=" , num2str(N) , "$)" );
legend3 = "Without KPM";
legend({legend1,legend2,legend3},'Interpreter','latex');
    

xlim([-1,+1]);
xlabel('$x$','Interpreter','latex');
xticks([ -1.00 , -0.75 , -0.50 , -0.25 , 0.00 , +0.25 , +0.50 , +0.75 , +1.00 ]);

ylim([-5,21]);
ylabel('$y$','Interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on');
