% MATLAB R2018a
clear;
clc;
close all;
addpath('../functions');


figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,1,1.0,1.0]);


Nsizes = [ 16 , 32 , 64 , 128 , 256 ];
for i = 1 : length(Nsizes)
    
    N = Nsizes(i);
    
    lambda = 1.5;
    mun = cos( pi/2 * [0:N-1] ) ./ ( 1 + ([0:N-1]==0) );
    mun = function_applyKernel( mun , 'Lorentz' , lambda );
    
    
    gamma = lambda/N;
    xPlotMesh1 = linspace(-0.5,+0.5,1000);
    yPlotMesh1 = (gamma/pi) ./ ( gamma^2 + xPlotMesh1.^2 );
   
    [xPlotMesh2,yPlotMesh2] = function_fftEvaluateChebSeries( mun , 1.0 , 100000 );
    
 
    subplot(1,length(Nsizes),i);
    plot( xPlotMesh1 , yPlotMesh1 , 'ro' , 'LineWidth' , 2.0 , 'MarkerSize' , 3 ); hold on;
    plot( xPlotMesh2 , yPlotMesh2 , 'b-' , 'LineWidth' , 1.0 , 'MarkerSize' , 3 ); hold on;
    
    legend1 = 'Lorentzian $\gamma=\lambda/N$';
    legend2 = 'KPM ($\lambda= ' + string(lambda) + ',N=' + string(N) + ')$';
    legend({legend1,legend2},'Interpreter','latex','FontSize',10);
    
    grid on; grid minor;
    xlim([-0.5,+0.5]);ylim([0,55]);
    xlabel('$x$','Interpreter','latex','FontSize',30);
    if( i == 1 ) ylabel('$f(x)$','Interpreter','latex','FontSize',30); end
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',20);

end
